import os
import random
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import json
import scipy.stats as st
import sys

import tkinter as tk
from tkinter import Label, ttk, filedialog

from GUI.tooltip import ToolTip
from locus import Locus
from locus_tracker import LocusTracker
from global_variables import variables

from timeit import default_timer as timer

# Set global variables
cheater_allele_counts = []
resistor_allele_counts = []
mean_ch_eff = []
mean_res_eff = []
wild_counts = []

locus_plot = []

mt1_ratios = []
mt2_ratios = []
mt3_ratios = []

# Create labels and entry widgets for each variable
variables_tooltips = {
    "n_runs": "Number of repeat simulation runs",
    "n_dev": "Total development cycles",
    "i_macrocyst": "In which development cycles should sexual cycles occur? Formatted as python list",
    "sex_cycle_interval": "Undergo sexual cycles at this interval - i_macrocyst must equal 0 if this value is higher than 0 (ie its either the list or\n the interval, not both)",
    "vg": "Vegetative growth generations per sexual cycle",
    "n": "Population size (rounds up to nearest even number)",
    "sl": "Slug size",
    "sp": "Spores per slug (fraction)",
    "m_ch": "Mutation rate to cheater allele",
    "m_re": "Mutation rate to resistance allele",
    "c_ch": "Cheater allele cost during vegetative growth",
    "c_res": "Resistance allele cost during vegetative growth",
    "c_ch_res": "Combined cheater + resistance cost - deviation from multiplicative (deviation from expectation)",
    "ch_germ": "Cheater germination rate penalty, ie 1 = cheaters never germinate, 0 = cheaters always germinate",
    "ch_start": "Starting cheater allele frequency",
    "res_start": "Starting resistance allele frequency",
    "ch_res_dist": "Cheater and resistance alelles initially distributed randomly (1) or without co-occurence (0)",
    "ch_eff_rc": "Multiplier on ch_eff in loci with both alleles(float between 0 and 1, 0 = cheater nullified to 1 = no effect)",
    "res_eff_rc": "Multiplier on res_eff in loci with both alleles(float between 0 and 1, 0 = resistance nullified to 1 = no effect)",
    "mc_count": "Number of macrocysts created per cycle, number\nof cells cannibilized per macrocyst is dependent on total pop size",
    "mc_germ_count": "Number of new cells generated per macrocysts, ie 1 macrocyst germinates into X cells",
    "recomb_chance": "Chance of recombining parental genotypes",
    "mt1_start": "Ratio of mating type 1 at start.\nmt1 + mt2 + m3 must equal 1",
    "mt2_start": "Ratio of mating type 2 at start.\nmt1 + mt2 + m3 must equal 1",
    "mt3_start": "Ratio of mating type 2 at start.\nmt1 + mt2 + m3 must equal 1",
    "output_filepath": "Filepath to save the run parameters and output to (.txt).\nLeave blank if you don't want to save the output",
    "gene_pairs": "Number of loci to simulate",
    "confidence_interval": "Desired confidence interval displayed when graphing multiple runs"
}

# Defining GUI functions and classes

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

def submit_form():
    for var_name, entry_widget in entry_widgets.items():
        # Handle number variables, split into floats and integers
        if var_name != "i_macrocyst" and var_name != "output_filepath":
            if var_name in list(["n_runs", "n_dev", "vg", "n", "sl", "ch_res_dist", "mc_count", "mc_germ_count", "gene_pairs", "ch_self_cheat", "sex_cycle_interval", "ch_eff_rc", "res_eff_rc"]):
                value = entry_widget.get()
                variables[var_name] = int(value)
            else:
                value = entry_widget.get()
                variables[var_name] = float(value)
        # Handles list parameters
        elif var_name == "i_macrocyst":
            value = entry_widget.get()
            split_list = value.split(',')
            for character in split_list:
                character.strip()
            new_split_list = list(map(int, split_list))
            print(new_split_list)
            variables[var_name] = new_split_list
        # Handles str parameters (this should probably check if its a valid file path but it doesn't)
        elif var_name == "output_filepath":
            value = entry_widget.get()
            variables[var_name] = str(value)
    progressbar = ttk.Progressbar(variable = progress)
    progressbar.grid(row=15, column=3, columnspan=3, pady=10)

    main()

# Simulation
class Cell:
    _id_counter = 0  # Class attribute to keep track of the last ID assigned

    def __init__(self, positionX, positionY, fitness, mating_type, spore_chance, loci_list):
        self.loci = loci_list
        self.positionX = positionX
        self.positionY = positionY
        self.fitness = fitness
        self.mating_type = mating_type
        self.spore_chance = spore_chance

        self.key = Cell._id_counter
        Cell._id_counter += 1
    
    def __str__(self):
        return f"{self.key}"
    
    def __repr__(self):
        return f"{self.key}"
    
    def __eq__(self, other):
            if isinstance(other, Cell):
                return self.key == other.key
            return False
    def clone(self):
        # This function is needed to ensure each time a cell is duplicated it is assigned a new spot in memory - without it duplicated cells
        # point to the same object and a change to one changes all
        new_loci_list = [Locus(locus.cheater_allele, locus.resistor_allele) for locus in self.loci]
        return Cell(self.positionX, self.positionY, self.fitness, self.mating_type, self.spore_chance, new_loci_list)
    
    def __hash__(self):
        return hash(self.key)

    def cheater_value(self):
        value = 0
        for locus in self.loci:
            value += locus.cheater_allele
        return value

    def wild_value(self):
        value = 0
        for locus in self.loci:
            if locus.cheater_allele == 0:
                value += 1
            if locus.resistor_allele == 0:
                value += 1
        return value

    def resistor_value(self):
        value = 0
        for locus in self.loci:
            value += locus.resistor_allele
        return value

    def ch_eff(self):

        # Effectiveness of a single loci
        single_loci_value = 1/len(self.loci)

        # Count how many loci have both alleles present
        both_alleles = 0
        for locus in self.loci:
            both_alleles += locus.both()

        # Overall ch effectiveness is equal to # loci with only ch allele times single loci effect, plus (# loci with both weighted by ch_eff_rc)
        eff = ((self.cheater_value()-both_alleles)*single_loci_value)+((both_alleles*variables["ch_eff_rc"])*single_loci_value)
        return eff

    def res_eff(self):
        # Effectiveness of a single loci
        single_loci_value = 1/len(self.loci)

        # Count how many loci have both alleles present
        both_alleles = 0
        for locus in self.loci:
            both_alleles += locus.both()

        # Overall ch effectiveness is equal to # loci with only ch allele times single loci effect, plus (# loci with both weighted by ch_eff_rc)
        eff = ((self.resistor_value()-both_alleles)*single_loci_value)+((both_alleles*variables["res_eff_rc"])*single_loci_value)
        return eff
    
    def exploitable(self, other):

        #If the cheater allele is invalidated by co-occurence
        if variables["ch_eff_rc"] == 0:
            for index in range(0, len(self.loci)):
                if self.loci[index].both() == 1 or self.loci[index].cheater_allele == 0: #Check if there is co-occurence (invalidates ch allele) or if there is no cheater allele
                    continue # If cant cheat at this gene pair, go to next pair
                elif variables["res_eff_rc"] == 0: # If resistance allele is invalidated by co-occurence
                    if other.loci[index].cheater_allele == 1 and other.loci[index].both() == 0: # First check if the other cell has an active cheater allele at this pair
                        continue # If it does, cant cheat go to next pair
                    elif other.loci[index].both() == 1 or other.loci[index].resistor_allele == 0: # Check if other cell has co-occurence or no resistor (same thing here)
                        return 1 # If it does, this cell can cheat the other. Return 1 and exit
                elif variables["res_eff_rc"] == 1: # If the resistance allele does not have epistasis
                    if other.loci[index].cheater_allele == 1 and other.loci[index].both() == 0: # Check if the toher cell has an active cheater allele
                        continue # This pair cant be cheated, go to next pair
                    elif other.loci[index].resistor_allele == 0: #If the other cell does not have a resistor allele, can cheat.
                        return 1   
                 
        elif variables["ch_eff_rc"] == 1: #If epsistasis does not occur for cheater allele
            for index in range(0, len(self.loci)):
                if self.loci[index].cheater_allele == 0: # If this pair does not have a cheater allele
                    continue # Go to next pair, this pair cant cheat
                elif variables["res_eff_rc"] == 0: #If 
                    if other.loci[index].cheater_allele == 1:
                        continue
                    elif other.loci[index].both() == 1 or other.loci[index].resistor_allele == 0:
                        return 1
                elif variables["res_eff_rc"] == 1:
                    if other.loci[index].cheater_allele == 1:
                        continue
                    elif other.loci[index].resistor_allele == 0:
                        return 1  
                
        return 0 

def load_parameters():
    root.update()
    filename = filedialog.askopenfilename()
    root.update()

def validate_entry(inp):
    if inp == "":
        return True
    try:
        float(inp)
    except:
        return False
    return True


def fifty_fifty(a, b):
    if np.random.randint(0,2) == 1:
        return a
    else:
        return b   

def mean_confidence_interval(data, confidence=variables["confidence_interval"]):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), st.sem(a)
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    return h

#Start up the population, returns a list of cell objects
def initialise_pop():
    global locus_plot

    pop_list = []
    mating_type_list= [1, 2, 3]
    mating_type_chances = [variables["mt1_start"], variables["mt2_start"], variables["mt3_start"]]

    # Inititate blank loci list for creating new Cell objects and loci tracker list
    blank_loci_list = []
    locus_tracker_list = []

    for x in range(0, variables["gene_pairs"]):
        name = "Locus"+str(x)
        blank_loci_list.append(Locus(0,0))
        locus_tracker_list.append(LocusTracker(name, 0, 0, 0))

    print("length of loci: " + str(len(blank_loci_list)))


    #Make a new list of length variables["n"] filled with cell objects
    for i in range(0, int(variables["n"])):

        chosen_mt = random.choices(mating_type_list, k=1, weights=mating_type_chances)[0]
        #Create the cell objects with default values and random mating type
        new_list = []
        for locus in blank_loci_list:
            new_list.append(locus.clone())
        pop_list.append(Cell(0, 0, 1, chosen_mt, 0, new_list.copy()))
    
    #Randomly select cells by index to add cheater/resistance alleles to
    if variables["ch_res_dist"] == 1:
        cheater_index_list = random.sample(list(range(variables["n"])), round(variables["ch_start"]*variables["n"]))
        resistor_index_list = random.sample(list(range(variables["n"])), round(variables["res_start"]*variables["n"]))
    else:
        cheater_index_list=[]
        resistor_index_list=[]
        ch_end = int((variables["ch_start"]*variables["n"]))
        last_index = 0
        for i in range(0, ch_end):
            cheater_index_list.append(i)
            last_index = i

        res_end = int(last_index + (variables["res_start"]*variables["n"]))
        for i in range(last_index, res_end):
            resistor_index_list.append(i)

    #Set the alleles at each locus, i.e. if cell was selected as cheater all loci have cheater allele, same for resistor
    check = len(cheater_index_list)
    for num in cheater_index_list:
        for locus in pop_list[num].loci:
            locus.cheater_allele = 1

    for num in resistor_index_list:
        for locus in pop_list[num].loci:
            locus.resistor_allele = 1

    
    global cheater_allele_counts
    global resistor_allele_counts
    global wild_counts
    global mean_ch_eff
    global mean_res_eff

    global mt1_ratios
    global mt2_ratios
    global mt3_ratios
    
    cheater_allele_count = 0
    wild_count = 0
    resistor_allele_count = 0
    total_ch_eff = 0
    total_res_eff = 0
    res_count = 0
    ch_count = 0

    mt1 = 0
    mt2 = 0
    mt3 = 0

    for cell in pop_list:
        cheater_allele_count += cell.cheater_value()
        resistor_allele_count += cell.resistor_value()
        wild_count += cell.wild_value()

        locus_num = 0
        for locus in cell.loci:
            if locus.both() == 0:
                locus_tracker_list[locus_num].ch_count += locus.cheater_allele
                locus_tracker_list[locus_num].res_count += locus.resistor_allele
            else:
                locus_tracker_list[locus_num].both_count += 1
            locus_num += 1

        if cell.cheater_value() > 0:
            total_ch_eff += cell.ch_eff()
            ch_count += 1
        if cell.resistor_value() > 0:
            total_res_eff += cell.res_eff()
            res_count += 1

        if cell.mating_type == 1:
            mt1 += 1
        elif cell.mating_type == 2:
            mt2 += 1
        elif cell.mating_type == 3:
            mt3 += 1

    locus_snapshot = []
    for tracker in locus_tracker_list:
        locus_snapshot.append(tracker.clone())
    locus_plot.append(locus_snapshot)

    total_open_slots = variables["gene_pairs"]*len(pop_list)
    
    cheater_allele_counts.append(cheater_allele_count/total_open_slots)
    resistor_allele_counts.append(resistor_allele_count/total_open_slots)
    wild_counts.append(wild_count/(total_open_slots*2))

    # Avoid divide by zero error
    if ch_count == 0:
        mean_ch_eff.append(0)
    else:
        mean_ch_eff.append(total_ch_eff/ch_count)
    if res_count == 0:
        mean_res_eff.append(0)
    else:
        mean_res_eff.append(total_res_eff/res_count)

    mt1_ratios.append(mt1/len(pop_list))
    mt2_ratios.append(mt2/len(pop_list))
    mt3_ratios.append(mt3/len(pop_list))

    return pop_list

#Vegatative growth cycle
def vegetative_growth(pop_list):
    pop_copy = pop_list
    kill_list =[]

    #Loop through the population
    for index in range(0,len(pop_copy)):
        #Check if spore is a cheater
        if pop_copy[index].cheater_value() > 0:
            #If it is, roll a number between 0-100, if the number is greater than the cheater chance to germinate (variables["ch_germ"]) add the index to the list of spores to be killed
            if (1-(variables["ch_germ"]*pop_copy[index].ch_eff())) < random.uniform(0, 1):
                kill_list.append(index)

    # Delete spores in kill_list from main pop
    for index in sorted(kill_list, reverse=True):
        del pop_copy[index]

    #Loop for number of vegetative growth cycles
    for j in range(0, variables["vg"]):
        new_pop=[]
        #Loop through the population
        for cell in pop_copy:
            # Set cell fitness to the average of cheater/resistor fitness, weighted by effectiveness of each
                cell.fitness = ((1-(variables["c_ch"] * cell.cheater_value()))) + (1-(variables["c_res"] * cell.resistor_value()))/2


        # Extract fitness values from Cell objects
        fitness_values = [cell.fitness for cell in pop_copy]

        # Sample with replacement based on fitness, ie fitness = chance to be selected
        selected_indices = random.choices(range(len(pop_copy)), k=variables["n"], weights=fitness_values)

        # Create the next population generation based on replicating individuals
        new_pop = [pop_copy[i].clone() for i in selected_indices]
        pop_copy = new_pop

        # Add new mutations - invert current allele if the allele rolls under the mutation chance
        # Can remove exisiting alleles
        for cell in pop_copy:
            for locus in cell.loci:
                if (random.uniform(0, 1) <= variables["m_ch"]):
                    if locus.cheater_allele == 1:
                        locus.cheater_allele = 0
                    else:
                        locus.cheater_allele = 1

                if (random.uniform(0, 1) <= variables["m_re"]):
                    if locus.resistor_allele == 1:
                        locus.resistor_allele = 0
                    else:
                        locus.resistor_allele = 1

    
    #Return the new population of cells
    return pop_copy

def developmental_cycle(dpop_list):
    global locus_plot

    #Empty list to store all cells that turn into spores from respective slugs
    fruiting_body = []
    total_stalk = 0

    #Loop through each formable slug based on the size of the pop and the size of the slugs - leftover cells are ignored
    for i in range(0,math.floor(variables["n"]/variables["sl"])):
        #Make an empty slug
        slug = []

        #Select 'variables["sl"]' number of cells from the pop
        slug_cell_indexes = random.sample(list(range(len(dpop_list))), variables["sl"])
        #Add selected cells to the slug
        for i in slug_cell_indexes:
            slug.append(dpop_list[i].clone())

        #Delete cells that formed slug from total pop
        for index in sorted(slug_cell_indexes, reverse=True):
            del dpop_list[index]


        ################################################
        # Determine spore/stalk cell fates
        # Stalk first determination
        ################################################

        # Initiate lists outside of loop
        stalk_index_list = []
        pre_stalk_list = []
        slug_list = []

        for i in range(0, len(slug)):
            slug_list.append(i)

        for i in range(0, len(slug)):
            if random.uniform(0, 1) < (1-variables["sp"]):
                pre_stalk_list.append(i)

        for num in pre_stalk_list:
            if slug[num].cheater_value() > 0:
                set_stalk_index = set(stalk_index_list)
                set_slug = set(slug_list)
                set_prestalk = set(pre_stalk_list) 
                result_set = set_stalk_index.symmetric_difference(set_slug)
                final_result_set = result_set.symmetric_difference(set_prestalk)
                result_list = list(final_result_set)

                other_cell_index = result_list[np.random.randint(0, len(result_list))]

                if slug[num].exploitable(slug[other_cell_index]) == 1:
                    stalk_index_list.append(other_cell_index)
                else:
                    stalk_index_list.append(num)
            else:
                stalk_index_list.append(num)

        # Remove stalk cells from slug (iterate through indexes backwards)
        total_stalk = total_stalk + len(stalk_index_list)
        for index in sorted(stalk_index_list, reverse=True):
            del slug[index]

        fruiting_body.extend(slug)

    locus_tracker_list = []
    loci_track_count = 0
    for x in range(0, variables["gene_pairs"]):
        name = "Locus"+str(loci_track_count)
        locus_tracker_list.append(LocusTracker(name, 0, 0, 0))
        loci_track_count += 1
    
    cheater_allele_count = 0
    wild_count = 0
    resistor_allele_count = 0
    total_res_eff = 0
    total_ch_eff = 0
    ch_count = 0
    res_count = 0


    mt1 = 0
    mt2 = 0
    mt3 = 0

    for cell in fruiting_body:
        cheater_allele_count = cheater_allele_count + cell.cheater_value()
        resistor_allele_count = resistor_allele_count + cell.resistor_value()
        wild_count = wild_count + cell.wild_value()

        locus_num = 0
        for locus in cell.loci:
            if locus.both() == 0:
                locus_tracker_list[locus_num].ch_count += locus.cheater_allele
                locus_tracker_list[locus_num].res_count += locus.resistor_allele
            else:
                locus_tracker_list[locus_num].both_count += 1
            locus_num += 1

        if cell.cheater_value() > 0:
            total_ch_eff += cell.ch_eff()
            ch_count += 1
        if cell.resistor_value() > 0:
            total_res_eff += cell.res_eff()
            res_count += 1
        
        if cell.mating_type == 1:
            mt1 += 1
        elif cell.mating_type == 2:
            mt2 += 1
        elif cell.mating_type == 3:
            mt3 += 1

    total_open_slots = variables["gene_pairs"]*len(fruiting_body)
    
    print("Len: " +str(len(fruiting_body))  + ", cheater: " + str(cheater_allele_count) + ", resistor: " + str(resistor_allele_count) + ", wild: " + str(wild_count) + " , total stalks: " + str(total_stalk))

    global cheater_allele_counts
    global resistor_allele_counts
    global wild_counts

    global mean_ch_eff
    global mean_res_eff

    global mt1_ratios
    global mt2_ratios
    global mt3_ratios

    cheater_allele_counts.append(cheater_allele_count/total_open_slots)
    resistor_allele_counts.append(resistor_allele_count/total_open_slots)
    wild_counts.append(wild_count/(total_open_slots*2))

    locus_snapshot = []
    for tracker in locus_tracker_list:
        locus_snapshot.append(tracker.clone())
    locus_plot.append(locus_snapshot)

    if ch_count == 0:
        mean_ch_eff.append(0)
    else:
        mean_ch_eff.append(total_ch_eff/ch_count)
    if res_count == 0:
        mean_res_eff.append(0)
    else:
        mean_res_eff.append(total_res_eff/res_count)


    mt1_ratios.append(mt1/len(fruiting_body))
    mt2_ratios.append(mt2/len(fruiting_body))
    mt3_ratios.append(mt3/len(fruiting_body))

    if("--param" not in  sys.argv):
        global progress
        progress.set(progress.get() + (100/variables["n_dev"]))
        root.update()

    return fruiting_body

def sexual_cycle(sexual_pop):
    sex_pop = sexual_pop
    mc_pop = []
    first_selection_indexes = random.sample(list(range(len(sex_pop))), variables["mc_count"])
    return_pop = []

    # Add selected indexes to new list (these will be used as the founders of each macrocyst) and remove them from the old list
    for num in first_selection_indexes:
        mc_pop.append(sex_pop[num])

    sex_pop = [item for item in sex_pop if item not in mc_pop]

    # Create empty loci list for progeny
    blank_list = []
    for i in range(0, variables["gene_pairs"]):
        blank_list.append(Locus(0,0)) 

    # Loop through every founder cell, pick a partner, generate progeny and populate the return list
    for cell in mc_pop:
        temp_list = []
        # Make a list of all available partners in the original population, pick one at random and remove it from the population
        for pop_cell in sex_pop:
            if cell.mating_type != pop_cell.mating_type:
                temp_list.append(pop_cell)
        partner_cell = random.sample(temp_list, 1)[0]
        sex_pop.remove(partner_cell)

        # Initiate empty progeny cell
        progeny_cell = Cell(0, 0, 1, 0, 0, blank_list.copy())

        # If both mating cells share the same cheating/resistance genes, so does the progeny
        if cell.loci == partner_cell.loci:
            progeny_cell.loci = cell.loci.copy()
            progeny_cell.mating_type = fifty_fifty(cell.mating_type, partner_cell.mating_type)

        # Roll a number to decide if the progeny will undergo recombination
        # No recombination
        elif random.uniform(0, 1) < (1-variables["recomb_chance"]):
            #If no recombination, progeny will have one of the parental genotypes at this locus (equal chance)
            progeny_cell.loci = fifty_fifty(cell.loci, partner_cell.loci)

            # Give the progeny cell the parental mating type depending on parent chosen for locus
            if progeny_cell.loci == cell.loci:
                progeny_cell.mating_type = cell.mating_type
            else:
                progeny_cell.mating_type = partner_cell.mating_type

        # Recombination
        else:
            count = 0
            for locus in progeny_cell.loci:
                locus.cheater_allele = fifty_fifty(cell.loci[count].cheater_allele, partner_cell.loci[count].cheater_allele)
                locus.resistor_allele = fifty_fifty(cell.loci[count].resistor_allele, partner_cell.loci[count].resistor_allele)
                count += 1
            progeny_cell.mating_type = fifty_fifty(cell.mating_type, partner_cell.mating_type)
        
        #Set fitness based on loci
        progeny_cell.fitness = ((1-(variables["c_ch"] * progeny_cell.ch_eff())) + (1-(variables["c_res"] * progeny_cell.res_eff())))/2

        #Add variables["mc_germ_count"] number of progeny to the returned population
        for i in range(0, variables["mc_germ_count"]):
            return_pop.append(progeny_cell.clone())

    print("Completed sexual cycle")
    return return_pop

def main():
    # Get global variables
    global cheater_allele_counts
    global resistor_allele_counts
    global wild_counts
    global mean_ch_eff
    global mean_res_eff
    global mt1_ratios
    global mt2_ratios
    global mt3_ratios
    global locus_plot

    # Initiate graphing lists
    total_ch = []
    mean_ch = []
    ch_ci = []

    total_res = []
    mean_res = []
    res_ci = []

    total_wild = []
    mean_wild = []
    wild_ci = []

    total_mt1 = []
    mean_mt1 = []
    mt1_ci = []

    total_mt2 = []
    mean_mt2 = []
    mt2_ci = []

    total_mt3 = []
    mean_mt3 = []
    mt3_ci = []

    graph_ch_eff = []
    graph_mean_ch_eff = []
    ch_eff_ci = []

    graph_res_eff = []
    graph_mean_res_eff = []
    res_eff_ci = []

    # Create tracker dictionary
    tracker_dict = {}
    average_tracker_dict = {}
    ci_tracker_dict = {}

    for i in range(0, variables["gene_pairs"]):
        temp_str = "gene_pair" + str(i)
        tmp1 = []
        tmp2 = []
        tmp3 = []
        tracker_dict[temp_str] = [tmp1.copy(), tmp2.copy(), tmp3.copy()]
        average_tracker_dict[temp_str] = [tmp1.copy(), tmp2.copy(), tmp3.copy()]
        ci_tracker_dict[temp_str] = [tmp1.copy(), tmp2.copy(), tmp3.copy()]

    # Initiate list of which rounds sexual cycle was performed
    sex_cycle_list = []

###################################################
    for j in range(0, variables["n_runs"]):
        # Initiliase and run round 0
        pop = initialise_pop()

        # Do the rest of the cycles
        for i in range(1, variables["n_dev"]+1):

            # If sexual cycles set by interval, check if current round is divisible by interval with no remainder - if so, do a sexual cycle
            if variables["sex_cycle_interval"] != 0 and (variables["i_macrocyst"][0] == 0):
                if i % variables["sex_cycle_interval"] == 0:
                    pop = sexual_cycle(pop)
                    pop = vegetative_growth(pop)
                    pop = developmental_cycle(pop)
                    sex_cycle_list.append(i)
                else:
                    pop = vegetative_growth(pop)
                    pop = developmental_cycle(pop)
                print("Dev round " + str(i) +" of " + str(variables["n_dev"])+ " completed.")
            
            # If sexual cycles set by list, check if current round is in the list and do a sexual cycle if it is
            else:
                if i in variables["i_macrocyst"]:
                    pop = sexual_cycle(pop)
                    pop = vegetative_growth(pop)
                    pop = developmental_cycle(pop)

                else:
                    pop = vegetative_growth(pop)
                    pop = developmental_cycle(pop)
                print("Dev round " + str(i) +" of " + str(variables["n_dev"])+ " completed.")

        # Locus tracking
        for num in range(0, len(locus_plot)):
            if j ==0:
                for i in range(0, variables["gene_pairs"]):
                    temp_str = "gene_pair" + str(i)
                    tracker_dict[temp_str][0].append([locus_plot[num][i].ch_count])
                    tracker_dict[temp_str][1].append([locus_plot[num][i].res_count]) 
                    tracker_dict[temp_str][2].append([locus_plot[num][i].both_count])
            else:
                for i in range(0, variables["gene_pairs"]):
                    temp_str = "gene_pair" + str(i)
                    tracker_dict[temp_str][0][num].append(locus_plot[num][i].ch_count)
                    tracker_dict[temp_str][1][num].append(locus_plot[num][i].res_count) 
                    tracker_dict[temp_str][2][num].append(locus_plot[num][i].both_count)

        # Begin building means/confidence interval lists for graphing
        for num in range(0, len(cheater_allele_counts)):
            if j == 0:
                # Genotype variables
                temp_list = [cheater_allele_counts[num]]
                total_ch.append(temp_list)
                temp_list = [resistor_allele_counts[num]]
                total_res.append(temp_list)
                temp_list = [wild_counts[num]]
                total_wild.append(temp_list)

                #Mating type variables
                temp_list = [mt1_ratios[num]]
                total_mt1.append(temp_list)
                temp_list = [mt2_ratios[num]]
                total_mt2.append(temp_list)              
                temp_list = [mt3_ratios[num]]
                total_mt3.append(temp_list)

                # Efficiency variables
                temp_list = [mean_ch_eff[num]]
                graph_ch_eff.append(temp_list)
                temp_list = [mean_res_eff[num]]
                graph_res_eff.append(temp_list)

            else:
                # Populate genotype graph lists
                total_ch[num].append(cheater_allele_counts[num])
                total_res[num].append(resistor_allele_counts[num])
                total_wild[num].append(wild_counts[num])

                # Populate mating-type graphing lists
                total_mt1[num].append(mt1_ratios[num])
                total_mt2[num].append(mt2_ratios[num])
                total_mt3[num].append(mt3_ratios[num])

                #Populate efficiency graphing lists
                graph_ch_eff[num].append(mean_ch_eff[num])
                graph_res_eff[num].append(mean_res_eff[num])

        # Reset trackers
        cheater_allele_counts = []
        resistor_allele_counts = []
        wild_counts = []
        mean_ch_eff = []
        mean_res_eff = []
        mt1_ratios = []
        mt2_ratios = []
        mt3_ratios = []
        locus_plot = []

        print("------------------")
        print("Run " + str(j+1) + " of " + str(variables["n_runs"])+" completed")
        print("------------------")

        if("--param" not in  sys.argv):
            # Reset loading bar
            global progress
            progress.set(0)
            root.update()

    # After all runs are finished, calculate means and confidence intervals
    x_axis_values = list(range(0, variables["n_dev"]+1))

    # Calculate the average # of cheater only, resistor only, and both at each gene pair for each time point and save to dict - this dict
    # will be saved to json and mapped with grapher.py
    for i in range(0, len(x_axis_values)):
        for z in range(0, variables["gene_pairs"]):
            temp_str = "gene_pair" + str(z)
            average_tracker_dict[temp_str][0].append(np.mean(tracker_dict[temp_str][0][i]))
            average_tracker_dict[temp_str][1].append(np.mean(tracker_dict[temp_str][1][i]))
            average_tracker_dict[temp_str][2].append(np.mean(tracker_dict[temp_str][2][i]))

            ci_tracker_dict[temp_str][0].append(mean_confidence_interval(tracker_dict[temp_str][0][i]))
            ci_tracker_dict[temp_str][1].append(mean_confidence_interval(tracker_dict[temp_str][1][i]))
            ci_tracker_dict[temp_str][2].append(mean_confidence_interval(tracker_dict[temp_str][2][i]))

    for i in range(0, len(x_axis_values)):
        # Genotype
        mean_ch.append(np.mean(total_ch[i]))
        mean_res.append(np.mean(total_res[i]))
        mean_wild.append(np.mean(total_wild[i]))

        ch_ci.append(mean_confidence_interval(total_ch[i]))
        res_ci.append(mean_confidence_interval(total_res[i]))
        wild_ci.append(mean_confidence_interval(total_wild[i]))

        # Mating-type
        mean_mt1.append(np.mean(total_mt1[i]))
        mean_mt2.append(np.mean(total_mt2[i]))
        mean_mt3.append(np.mean(total_mt3[i]))

        mt1_ci.append(mean_confidence_interval(total_mt1[i]))
        mt2_ci.append(mean_confidence_interval(total_mt2[i]))
        mt3_ci.append(mean_confidence_interval(total_mt3[i]))

        # Efficiency
        graph_mean_ch_eff.append(np.mean(graph_ch_eff[i]))
        graph_mean_res_eff.append(np.mean(graph_res_eff[i]))

        ch_eff_ci.append(mean_confidence_interval(graph_ch_eff[i]))
        res_eff_ci.append(mean_confidence_interval(graph_res_eff[i]))

    #Dump to output
    if variables["output_filepath"] != "":
        out = variables["output_filepath"] + ".json"

        value_dict = {
            "parameters": variables,
            "sex_cycle_list": sex_cycle_list,
            "x_axis_values": x_axis_values, 
            "mean_ch": mean_ch,
            "ch_ci": ch_ci,
            "mean_res": mean_res,
            "res_ci": res_ci,
            "mean_wild": mean_wild,
            "wild_ci": wild_ci,
            "mean_mt1": mean_mt1,
            "mt1_ci": mt1_ci,
            "mean_mt2": mean_mt2, 
            "mt2_ci": mt2_ci, 
            "mean_mt3": mean_mt3,
            "mt3_ci": mt3_ci,
            "graph_mean_ch_eff": graph_mean_ch_eff,
            "ch_eff_ci": ch_eff_ci,
            "graph_mean_res_eff": graph_mean_res_eff,
            "res_eff_ci": res_eff_ci,
            "average_tracker_dict": average_tracker_dict,
            "ci_tracker_dict": ci_tracker_dict
        }
        with open(out, "w") as outfile: 
            json.dump(value_dict, outfile)


    if("--param" not in  sys.argv):
        # Plot graph
        fig, ax = plt.subplots(2) # change to 3 for efficacy
        fig.set_size_inches(12, 8)

        # Make patches for legends
        ch_patch = mpatches.Patch(color='red', label='Cheater alleles')
        re_patch = mpatches.Patch(color='blue', label='Resistor alleles')
        wild_patch = mpatches.Patch(color='black', label='Wild-type alleles')

        mt1_patch = mpatches.Patch(color='red', label='Type 1')
        mt2_patch = mpatches.Patch(color='blue', label='Type 2')
        mt3_patch = mpatches.Patch(color='green', label='Type 3')

        # Plot genotypes
        ax[0].set_ylabel("Ratio")
        ax[0].set_title("Allele frequency over time")
        ax[0].set_ylim([0, 1])

        ax[0].plot(x_axis_values, mean_ch, color = "red")
        ax[0].fill_between(x_axis_values, (np.array(mean_ch)-np.array(ch_ci)), (np.array(mean_ch)+np.array(ch_ci)), color='red', alpha=.1)
        ax[0].plot(x_axis_values, mean_res, color = "blue")
        ax[0].fill_between(x_axis_values, (np.array(mean_res)-np.array(res_ci)), (np.array(mean_res)+np.array(res_ci)), color='blue', alpha=.1)
        ax[0].plot(x_axis_values, mean_wild, color = "black")
        ax[0].fill_between(x_axis_values, (np.array(mean_wild)-np.array(wild_ci)), (np.array(mean_wild)+np.array(wild_ci)), color='black', alpha=.1)
        ax[0].grid()
        ax[0].legend(handles=[ch_patch, re_patch, wild_patch])

        # Plot mating types
        ax[1].set_ylabel("Ratio")


        ax[1].set_ylim([0, 1])

        ax[1].plot(x_axis_values, mean_mt1, color = "red")
        ax[1].plot(x_axis_values, mean_mt2, color = "blue")
        ax[1].plot(x_axis_values, mean_mt3, color = "green")
        ax[1].fill_between(x_axis_values, (np.array(mean_mt1)-np.array(mt1_ci)), (np.array(mean_mt1)+np.array(mt1_ci)), color='red', alpha=.1)
        ax[1].fill_between(x_axis_values, (np.array(mean_mt2)-np.array(mt2_ci)), (np.array(mean_mt2)+np.array(mt2_ci)), color='blue', alpha=.1)
        ax[1].fill_between(x_axis_values, (np.array(mean_mt3)-np.array(mt3_ci)), (np.array(mean_mt3)+np.array(mt3_ci)), color='green', alpha=.1)
        ax[1].grid()
        ax[1].legend(handles=[mt1_patch, mt2_patch, mt3_patch])

        # # Plot efficiency
        ax[0].set_xlabel("Development cycles")
        ax[1].set_xlabel("Development cycles")

        # Plot lines indicating sexual cycles
        if variables["i_macrocyst"][0] != 0:
            for i in variables["i_macrocyst"]:
                ax[0].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
                ax[1].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)

        else:
            for i in sex_cycle_list:
                ax[0].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)
                ax[1].axvline(x = i, color = "b", linestyle = "dashed", alpha = 0.05)

        plt.show()

###################
# Generating the GUI
####################

# Create the main window
if __name__ == "__main__":
    param = ""
    if("--param" in  sys.argv):
        param = sys.argv[sys.argv.index("--param") + 1]

        f = open(param)        
        variables = json.load(f)
        f.close()
        main()

    else:

        root = tk.Tk()
        root.title("Dictyostelium Cheater Simulation")

        progress = tk.IntVar()

        entry_widgets = {}
        row_num = 2
        column_num = 0

        # Create entry boxes for all items in variables
        for var_name, var_value in variables.items():
            label = tk.Label(root, text=var_name)
            label.grid(row=row_num, column=column_num, sticky="e")
            
            if var_name != "i_macrocyst" and var_name != "output_filepath":
                entry = tk.Entry(root, validate='key', vcmd=(root.register(validate_entry), '%P'))
                entry.insert(0, str(var_value))
                entry.grid(row=row_num, column=column_num+1)
                
            else:
                entry = tk.Entry(root)
                entry.insert(0, str(var_value))
                entry.grid(row=row_num, column=column_num+1)       

            # Give the entry boxes tooltips based on vairables_tooltips dictionary
            CreateToolTip(entry, text = variables_tooltips[var_name])

            entry_widgets[var_name] = entry
            row_num += 1
            if row_num + 2 >= 11:
                column_num = column_num + 2
                row_num = 2

        # Create a button to submit the form
        submit_button = tk.Button(root, text="Run simulation", command=submit_form)
        submit_button.grid(row=14, column=3, columnspan=3, pady=10)

        # Create a button to load parameters from txt file
        load_button = tk.Button(root, text="Import parameters", command=load_parameters)
        load_button.grid(row=14, column=1, columnspan=3, pady=10)
        load_button["state"]="disabled"

        # Create a button to run a suite of tests to ensure simulation works
        test_button = tk.Button(root, text="Run test suite", command=submit_form)
        test_button.grid(row=14, column=5, columnspan=3, pady=10)
        test_button["state"]="disabled"

        # Create the title label
        parameters_label = Label(root, text="Parameters")
        parameters_label.grid(row=0, column=3, columnspan=3, pady=10)

        print(f'Proccess ID: {os.getpid()}')

        # Start the Tkinter event loop
        root.mainloop()
