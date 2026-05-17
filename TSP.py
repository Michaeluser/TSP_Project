import math
import os
import random
import numpy
import pygame
import tkinter
from PIL import Image, ImageTk

#GENERAL_PARAMETERS
SOL_METHOD = "SA"
pind_amount = random.randint(20, 40) #number of genes in every permutation

MS_X = 700
MS_Y = 700


#VIZUALIZATION_PARAMETERS
BORDER = 50
HEIGHT = MS_Y + 2 * BORDER
WIDTH = MS_X + 2 * BORDER

#VISUALIZATION
def display_image(frame_pos, label, frame_list):
    img = Image.open(frame_list[int(frame_pos) - 1])
    tk_img = ImageTk.PhotoImage(img)
    label.config(image=tk_img)
    label.image = tk_img

def display_results(frame_list):
    root = tkinter.Tk()
    root.title(f"Results")
    root.geometry(f"{WIDTH}x{HEIGHT}")
    
    img = Image.open(frame_list[0])
    tk_img = ImageTk.PhotoImage(img)
    label = tkinter.Label(root, image = tk_img)
    
    slider = tkinter.Scale(root, from_= 1, to = len(frame_list), orient = "horizontal", command = lambda frame_pos: display_image(frame_pos, label, frame_list))
    

    label.pack(side = "top", expand = True)
    slider.pack(side = "bottom", fill = "x")

    root.mainloop()

def render(lr_permutations_list):
    pygame.init()
    screen = pygame.display.set_mode((WIDTH, HEIGHT))
    screen.fill((10, 10, 255))
    clock = pygame.time.Clock()
    
    for point in range(pind_amount):
        curr_ind = pr_permutations_list[point]
        pygame.draw.circle(screen, (255, 255, 255), (curr_ind[0] + BORDER, curr_ind[1] + BORDER), 4)
    
    if (SOL_METHOD == "GA"):
        pygame.display.set_caption(f"Visualization of Genetic Algorithm")
    else:
        pygame.display.set_caption(f"Vizualization of Simulated Annealing")
    
    pygame.display.flip()
    image_list = [-1] * len(lr_permutations_list)
    
    r_state = True
    perm_ordr = 0
    gen_len = len(lr_permutations_list)
    while r_state and perm_ordr < gen_len:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                r_state = False        

        screen.fill((10, 10, 255))
        permutation = lr_permutations_list[perm_ordr]
        
        for gene in range(pind_amount - 1):
            c_indiv = pr_permutations_list[permutation[gene]]
            n_indiv = pr_permutations_list[permutation[gene + 1]]
            pygame.draw.line(screen, (255, 255, 255), (c_indiv[0] + BORDER, c_indiv[1] + BORDER), (n_indiv[0] + BORDER, n_indiv[1] + BORDER), 2)
        
        pygame.draw.line(screen, (255, 255, 255), (pr_permutations_list[permutation[0]][0] + BORDER, pr_permutations_list[permutation[0]][1] + BORDER), (pr_permutations_list[permutation[pind_amount - 1]][0] + BORDER, pr_permutations_list[permutation[pind_amount - 1]][1] + BORDER), 2)
        
        clock.tick(60)
        pygame.display.flip()
        pygame.image.save(screen, f"frame_{perm_ordr + 1}.png")
        image_list[perm_ordr] = f"frame_{perm_ordr + 1}.png"
            
        perm_ordr += 1
        
    pygame.display.quit()
    display_results(image_list)


#GENERAL_UTILITY_FUNCTIONS
def save_data(data_list, file_name):
    global file  
    try:
        file = open(file_name, "a+")
        file.write("Length,Generation\n")
    except IOError:
        print("Couldn't read the file!")
    
    for datapoint in data_list:
        file.write(f"{datapoint[0]}, {datapoint[1]}\n")
    file.close()
    
    render(best_perm_chronology)
             


def gen_popul():
    pr_permutations_list = []
    i = 0
    while i < pind_amount:
        pr_permutations_list.append(tuple([random.randint(0, MS_X), random.randint(0, MS_Y)]))
        i += 1
    return pr_permutations_list


def perm_plength(permutation):
    plength = 0
    
    for k in range(pind_amount - 1):
        c_ind1 = pr_permutations_list[permutation[k]]
        c_ind2 = pr_permutations_list[permutation[k + 1]]
        
        c1 = (c_ind2[0] - c_ind1[0])
        c2 = (c_ind2[1] - c_ind1[1])
        
        plength += math.sqrt(c1**2 + c2**2)
    return plength



def define_segment(permutation, length):
    seg_pos = random.randint(0, 19 - length)
    seg_l = random.randint(1, length)
    seg = permutation[seg_pos:seg_pos + seg_l]
    [permutation.remove(el) for el in seg]
    return seg, seg_l, seg_pos, permutation


def mutate_permutation(permutation):
    perm_cpy = permutation.copy()
    
    mv = 7 #mutation variation amount
    probab_scale = [] # 1 - (displacement); 2 - (scramble); 3 - (inversion); 4 - (insertion); 5 - (adjacent swap); 6 - (random swap); 7 - (swap of 2 substrings)
    method_used = random.random()
    for p in range(mv):
        probab_scale.append(1/mv)
    probab_scale = numpy.cumsum(probab_scale)
    
    for prb in range(mv):
        if (probab_scale[prb] >= method_used):
            method_used = prb
            break
    match method_used:
        case 0:
            seg, seg_l, seg_pos, perm_cpy = define_segment(perm_cpy, 7)
            i = seg_l - 1
            while i > -1:
                perm_cpy.insert(seg_pos, seg[i])
                i -= 1
        case 1:
            seg, seg_l, seg_pos, perm_cpy = define_segment(perm_cpy, 7)
            random.shuffle(seg)
            i = seg_l - 1
            while i > -1:
                perm_cpy.insert(seg_pos, seg[i])
                i -= 1
        case 2:
            seg, seg_l, seg_pos, perm_cpy = define_segment(perm_cpy, 7)
            seg.reverse()
            i = seg_l - 1
            while i > -1:
                perm_cpy.insert(seg_pos, seg[i])
                i -= 1
            pass
        case 3:
            ipos = random.randint(0, pind_amount - 1)
            ind = perm_cpy[ipos]
            perm_cpy.remove(ind)
            perm_cpy.insert(random.randint(0, pind_amount - 1), ind)
        case 4:
            sswap_pos = random.randint(0, pind_amount - 2)
            ipop1 = perm_cpy[sswap_pos]
            ipop2 = perm_cpy[sswap_pos + 1]
            
            tmp = ipop1
            ipop1 = ipop2
            ipop2 = tmp
            del tmp
            
            perm_cpy[sswap_pos] = ipop1
            perm_cpy[sswap_pos + 1] = ipop2
        case 5:
            while True:
                s_pos1 = random.randint(0, pind_amount - 1)
                s_pos2 = random.randint(0, pind_amount - 1)
                
                if (s_pos1 != s_pos2):
                    ind1 = perm_cpy[s_pos1]
                    ind2 = perm_cpy[s_pos2]
                    tmp = ind2
                    
                    ind2 = ind1
                    ind1 = tmp
                    
                    perm_cpy[s_pos1] = ind1
                    perm_cpy[s_pos2] = ind2
                    break
        case 6:
            substring_length1 = random.randint(0, 8)
            substring_length2 = random.randint(0, 8)
            ss_pos1 = random.randint(0, int(pind_amount / 2 - substring_length1 - 1)) #swap starting position of 1st substring
            
            ss_pos2 = random.randint(ss_pos1 + substring_length1 + 1, pind_amount - 1) #swap starting position of 2nd substring
            while True:
                if not ((ss_pos2 + substring_length2) <= pind_amount - 1):
                    ss_pos2 = random.randint(ss_pos1 + substring_length1 + 1, pind_amount - 1)
                else:
                    break
            subs_1 = perm_cpy[ss_pos1:ss_pos1 + substring_length1 + 1]
            subs_2 = perm_cpy[ss_pos2:ss_pos2 + substring_length2 + 1]
            
            #ss_pos1 now represents relative position of substring 1, where subtring 2 will be inserted
            #ss_pos2 now represents length from the end of the list to the end of the substring 2 
            ss_pos2 = (pind_amount - 1) - (ss_pos2 + substring_length2) #
            comb_arrys = subs_1 + subs_2
            [perm_cpy.remove(el) for el in comb_arrys]
            
            [perm_cpy.insert(ss_pos1, subs_1[ind1]) for ind1 in range(len(subs_1) - 1, -1, -1)]
            [perm_cpy.insert(ss_pos2, subs_2[ind2]) for ind2 in range(len(subs_2))]
    
    return perm_cpy
    
    

#------------------------------GENETIC_ALGORITHM_SOLUTION_CODE------------------------------


#CONFIGURATION_VARIABLES_GA
GEN_NUM_GA = 700
PAR_SEL_METHOD = "RS" #Either SUS or RS
    
PS_PER_GENERATION_GA = 1000
WINDOW_SIZE = 20
MAX_STASIS_QUANTITY_GA = -1

if (PAR_SEL_METHOD == "RS"):
    MAX_STASIS_QUANTITY_GA = 200 #200
elif (PAR_SEL_METHOD == "SUS"):
    MAX_STASIS_QUANTITY_GA = 300 #300


#GENE_RECOMBINATION_AND_MUTATION
def recombine_permutations(perm1, perm2):
    dict_out = dict() 
    dict2 = dict()
    offspring = []
    for gpos in range(pind_amount):
        try:
            nb21 = perm1[gpos + 1]
            nb22 = perm2[gpos + 1]
        except IndexError:
            nb21 = perm1[0]
            nb22 = perm2[0]
        dict_out[perm1[gpos]] = {perm1[gpos - 1], nb21}
        dict2[perm2[gpos]] = {perm2[gpos - 1], nb22}
        
    for key1 in list(dict_out.keys()):
        for key2 in list(dict2.keys()):
            if (key1 == key2):
                dict_out[key1].update(dict2[key2])
    del dict2
    curr_gene = perm1[0]
    if (random.randint(0, 1) == 1):
        curr_gene = perm2[0]
    gene = 0
    while gene < pind_amount:
        for key in list(dict_out.keys()):
            dict_out[key].discard(curr_gene)
        offspring.append(curr_gene)
        
        if (len(dict_out[curr_gene])):
            neighbors_list = list(dict_out[curr_gene]) 
            stub_gene = neighbors_list[0]
            for curr_gene_neighbor in neighbors_list:
                if (len(dict_out[curr_gene_neighbor]) < len(dict_out[stub_gene])):
                    stub_gene = curr_gene_neighbor
            curr_gene = stub_gene
        else:
            for candidate in list(dict_out.keys()):
                if (candidate not in offspring):
                    curr_gene = candidate
                    break    
        gene += 1
        
    return offspring


#PARENT_SELECTION
def SUS_pselection(perm_plength_list, overall_fitness, parents_number, lr_perm_list): #lr_perm_list stands for logical representation permutation list
    probabilities_list = [-1] * PS_PER_GENERATION_GA
    for k in range(PS_PER_GENERATION_GA):
        probabilities_list[k] = (1/perm_plength_list[k])/overall_fitness
    probabilities_list = numpy.cumsum(probabilities_list)
    
    init_probability = 1 / random.randint(parents_number + 1,  2**20)
    parents = [-1] * parents_number
    
    for parent_pos in range(parents_number):
        parents[parent_pos] = init_probability + parent_pos/parents_number #here parents' list stores each parent's probability
    
    
    parent = 0
    while parent < parents_number:
        curr_parent_probability = parents[parent]
        for i in range(PS_PER_GENERATION_GA):
            if (probabilities_list[i] >= curr_parent_probability):
                parents[parent] = lr_perm_list[i] #here parents' list stores all parents
                break
        parent += 1
    return parents
    
def RS_pselection(perm_plength_list, parents_number, lr_perm_list): #lr_perm_list stands for logical representation permutation list
    probabilities_list = [-1] * PS_PER_GENERATION_GA
    parents = [-1] * parents_number
    perm_pll = sorted(perm_plength_list.copy())
    
    for i in range(PS_PER_GENERATION_GA):
        probabilities_list[i] = (2 * (PS_PER_GENERATION_GA - i))/(PS_PER_GENERATION_GA * (PS_PER_GENERATION_GA + 1))
    probabilities_list = numpy.cumsum(probabilities_list)
    init_probability = 1 / random.randint(parents_number + 1, 2**20)
    
    for parent_pos in range(parents_number):
        parents[parent_pos] = init_probability + parent_pos/parents_number #here parents' list contains each parent's probability
    
    
    parent = 0
    while parent < parents_number:
        curr_parent_probability = parents[parent]
        for prob_pos in range(PS_PER_GENERATION_GA):
            if (probabilities_list[prob_pos] >= curr_parent_probability):
                parents[parent] = lr_perm_list[perm_plength_list.index(perm_pll[prob_pos])]
                break
        parent += 1
    return parents
    

#ADDITIONAL_FUNCTIONS
def gen_permlist(perm_number):
    perm_list = []
    i = 0
    initial_perm_copy = initial_perm.copy()
    while i < perm_number:
        random.shuffle(initial_perm_copy)
        perm_list.append(initial_perm_copy.copy())
        i += 1
    return perm_list

def gen_fitness(gen_list):
    g_fitness = 0
    for perm in gen_list:
        g_fitness += (1/perm_plength(perm))
    return g_fitness
    



pr_permutations_list = gen_popul() #list of physical permutations i.e physical locations of each permutation
initial_perm = [ind for ind in range(pind_amount)] #starting permutation of the 1st population
lr_permutations_list = gen_permlist(PS_PER_GENERATION_GA) #list of logical permutations i.e logical displacement of cities in each permutation
perm_plength_list = [-1] * PS_PER_GENERATION_GA #list of path lengths of every permutation
data_list_GA = []
best_perm_chronology = []
fittest_perm = perm_plength(lr_permutations_list[0])

generation = 0
stasis_amount = 0
while SOL_METHOD == "GA":
    overall_fitness = gen_fitness(lr_permutations_list) #overall fitness of current generation
    parents = []
    offspring = []
    data_list_GA.append((fittest_perm, generation + 1))
    
    
    if (len(data_list_GA) % WINDOW_SIZE == 0):
        begin = int(len(data_list_GA) / WINDOW_SIZE) - 1
        for i in range(begin * WINDOW_SIZE, generation - 1):
            condition = ((data_list_GA[i][0] - data_list_GA[i + 1][0]) / data_list_GA[0][0]) < 1/1000
            if (condition):
                stasis_amount += 1
            else:
                stasis_amount = stasis_amount - 4
        if (stasis_amount >= MAX_STASIS_QUANTITY_GA or generation >= GEN_NUM_GA):
            if (PAR_SEL_METHOD == "SUS"):
                save_data(data_list_GA, "GA_succession_SUS.txt")
            elif (PAR_SEL_METHOD == "RS"):
                save_data(data_list_GA, "GA_succession_RS.txt")
            # render(best_perm_chronology)
            print("Change rate is too insignificant!\nHalting process!")
            os.abort()
    
    
    stub_var = 0
    fittest_index = 0
    while stub_var < PS_PER_GENERATION_GA:
        perm_plength_list[stub_var] = perm_plength(lr_permutations_list[stub_var]) #count physical lengths of every permutation in the current generation
        if (fittest_perm >= perm_plength_list[stub_var]):
            fittest_perm = perm_plength_list[stub_var]
            fittest_index = stub_var
        stub_var += 1
    best_perm_chronology.append(lr_permutations_list[fittest_index])
    
    
    if (PAR_SEL_METHOD == "SUS"):
        parents = SUS_pselection(perm_plength_list, overall_fitness, int(PS_PER_GENERATION_GA/10), lr_permutations_list)
    elif (PAR_SEL_METHOD == "RS"):
        parents = RS_pselection(perm_plength_list, int(PS_PER_GENERATION_GA/10), lr_permutations_list)
    
    stub_var = 0
    pp_amnt = PS_PER_GENERATION_GA/10 #first parent amount for pairs
    while stub_var < int(pp_amnt):
        par1 = parents[stub_var]
        par2_c = 0 #represents 2nd parent's ordinal number
        
        while par2_c < 10:
            par2_ind = random.randint(0, int(pp_amnt - 1))
            if (par2_ind == stub_var):
                continue
            mutation_prob = random.randint(0, 1000)
            recombination_prob = random.randint(0, 1000)
            
            
            if (650 <= recombination_prob <= 1000):
                offspring_permutation = recombine_permutations(par1, parents[par2_ind])
            else:
                par2 = parents[par2_ind]
                if (perm_plength(par1) > perm_plength(par2)):
                    offspring_permutation = par2
                else:
                    offspring_permutation = par1
                    
            
            if (1 <= mutation_prob <= 100):
                offspring_permutation = mutate_permutation(offspring_permutation)
            
            offspring.append(offspring_permutation)
            par2_c += 1
        stub_var += 1
    
    lr_permutations_list = offspring.copy()
    generation += 1


#------------------------------SIMULATED_ANNEALING_SOLUTION_CODE------------------------------

#CONFIGURATION_VARIABLES_SA
GEN_NUM_SA = 15000
TEMPERATURE = 80 #starting temperature 
MAX_STASIS_QUANTITY_SA = 1000 #maximum number of unchanged generations
DECLINE_RATE = 0.98 #rate at which temperature should change


best_perm = initial_perm.copy() #current best permutation
best_plength = perm_plength(best_perm) #current best permutation's length
data_list_SA = []

sel_prob = 0
curr_stasis = 0
generation = 0
curr_perm = list()
curr_plength = 0

while SOL_METHOD == "SA":
    if (curr_stasis > MAX_STASIS_QUANTITY_SA or generation > GEN_NUM_SA):
        save_data(data_list_SA, "SA_succession.txt")
        # render(best_perm_chronology)
        break
    
    
    curr_perm = mutate_permutation(best_perm)
    curr_plength = perm_plength(curr_perm)
    delta = curr_plength - best_plength
    exponent = -delta / TEMPERATURE    
    
    if (delta < 0):
        curr_stasis -= 1
        best_plength = curr_plength
        best_perm = curr_perm
        best_perm_chronology.append(best_perm)
        data_list_SA.append((best_plength, generation + 1))
    else:
        if (exponent < -700):
            sel_prob = 0.0
        else:
            sel_prob = math.exp(exponent)
            
        prob = random.random()
        if (prob <= sel_prob):
            curr_stasis -= 1
            best_perm = curr_perm
            best_plength = curr_plength
            best_perm_chronology.append(best_perm)
            data_list_SA.append((best_plength, generation + 1))
        else:
            curr_stasis += 1
            
    
            
    generation += 1
    TEMPERATURE *= DECLINE_RATE
