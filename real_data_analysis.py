# -*- coding: utf-8 -*-
"""

@author: Tim Diedrich
"""

import scipy.stats as stats
import igraph
from igraph import *
from matplotlib import cm
import math
from statistics import mean, variance


dataset_list = ["GSE3467", "GSE3678", "GSE4107", "GSE4183", "GSE8762",
                   "GSE14924_CD4", "GSE16515", "GSE16759", "GSE19188", "GSE19420", "GSE19728",
                   "GSE20153", "GSE22780", "GSE24739_G0", "GSE24739_G1",
                   "GSE32676", "GSE38666_epithelia", "GSE38666_stroma"]

basepath = "data_archive/"


#initial values for dataset GSE24739_G0 analysis
FILEPATH_ASSIGNMENT_MATRIX = basepath + "GSE24739_G0/assignment_matrix_filtered_min0_max100000.txt"
FILEPATH_KEGG = "data_archive/kegg.gmt"
FILEPATH_MALACARD_KEGG = basepath + "GSE24739_G0/malacard_kegg.csv"
FILEPATH_QVALUES_GENEID = basepath + "GSE24739_G0/GSE24739_G0_qvalues.csv"
FILEPATH_RES_PADOG = basepath + "GSE24739_G0/res_padog.csv"

FILEPATH_OUTPUT_JOANAPY = basepath + "GSE24739_G0/output_joanapy_single_species.csv"
FILEPATH_QVALUES_JOANA = basepath + "GSE24739_G0/qvalues_first_species.txt"
FILEPATH_TERMS = basepath + "GSE24739_G0/terms.txt"

OUTPUT_PATH = "plots/"

THRESHOLD = .1 #threshold for qvalues to be considered significant in independent-gene approach like Fisher's exact test



def get_rgb_color(percentage:float)->tuple:
    """
    helper function to return a gradient of green color for node colors

    Parameters
    ----------
    percentage : float
        percentage to select green gradient by

    Returns
    -------
    tuple
        rgb color gradient tuple to color graph nodes.

    """
    g = cm.get_cmap('Greens')
    if percentage == 0.0:
        return g(0)
    else:
        return g(percentage*0.75 + 0.3)



def fileLoader(path:str)->list:
    """
    helper function to load a single file content and return it

    Parameters
    ----------
    path : str
        path to file to load.

    Returns
    -------
    content : list
        file content.

    """
    with open(path) as f:
        content = f.readlines()       
    return content


def loadAllFiles()->tuple:
    """
    loads all files needed for analysis of a single dataset and modifies their content to be analyzable

    Returns
    -------
    asMatrix : list
        gene - gene set assignment matrix.
    kegg : list
        kegg gene - gene set assignment matrix.
    malKegg : list
        malacards relevance score according to the kegg dataset.
    qValuesGeneId : list
        q-values with gene ids.
    joanaOutput : list
        gene sets marked active by JOANA.
    joanaQValues : list
        q-values used to feed JOANA.
    terms : list
        list of all terms, i.e., gene sets.
    res_padog : list
        gene sets marked active by PADOG.

    """
    asMatrix = fileLoader(FILEPATH_ASSIGNMENT_MATRIX)
    asMatrix = [t.strip().split(",") for t in asMatrix]
    asMatrix = [[int(value) for value in t] for t in asMatrix]
    
    kegg = fileLoader(FILEPATH_KEGG)
    kegg = [t[:-1].split("\t")[2:] for t in kegg] #eliminate \n and descriptor
    
    malKegg = fileLoader(FILEPATH_MALACARD_KEGG)
    malKegg = [t.strip().split(",") for t in malKegg[1:]]
    
    qValuesGeneId = fileLoader(FILEPATH_QVALUES_GENEID)[1:] #header dropped
    
    joanaOutput = fileLoader(FILEPATH_OUTPUT_JOANAPY) 
    #joanaOutput = [t.strip().split(",") if t.count(",") == 1 else [",".join(t.strip().split(",")[:2])]+[t.strip().split(",")[2]] for t in joanaOutput[1:]] #header droppeda and split values
    joanaOutput = [t.strip().split(",") if t.count(",") == 1 else t.strip()[1:].split("\",") for t in joanaOutput[1:]] #header droppeda and split values
    
    joanaQValues = fileLoader(FILEPATH_QVALUES_JOANA)
    joanaQValues = [float(t.strip()) for t in joanaQValues]
    
    terms = fileLoader(FILEPATH_TERMS)
    terms = [t[:-1] if t.count("\"") == 0 else t[1:-2] for t in terms]
    
    res_padog = fileLoader(FILEPATH_RES_PADOG)
    res_padog = [t.strip().split(",") if t.count(",") <= 4 else [",".join(t.strip().split(",")[:2])]+t.strip().split(",")[2:] for t in res_padog[1:]]
    res_padog.sort()
    #check if padog values are complete
    #checkPadogComplete(res_padog)
    
    return asMatrix, kegg, malKegg, qValuesGeneId, joanaOutput, joanaQValues, terms, res_padog


#unused function, to asses completeness of PADOG data in datasets
def checkPadogComplete(padog):
    
    kegg = fileLoader(FILEPATH_KEGG)
    kegg = [t[:-1].split("\t") for t in kegg]
    for i in range(len(kegg)):
        
        if padog[i][0][0] == "\"":
            if kegg[i][0][:8] != padog[i][0][1:9]:
                padog.insert(i, ["\""+kegg[i][0][:8]+"\"", "dummy", "dummy", "dummy", 1])
                
        else:
            if kegg[i][0][:8] != padog[i][0][:8]:
                padog.insert(i, [kegg[i][0][:8], "dummy", "dummy", "dummy", 1])

    

def keggOverlap(kegg:list)->list:
    """
    creates an overlap matrix of how many genes overlap between the kegg gene sets

    Parameters
    ----------
    kegg : list
        kegg gene - gene set assignment matrix.

    Returns
    -------
    overlapMatrix : list
        matrix of number of overlapping genes.

    """
    overlapMatrix = []
    for t1 in kegg:   
        numOL = []
        for t2 in kegg:
            numGenes = 0
            for entry in t1:                
                if entry in t2:
                    numGenes+=1            
            numOL.append(numGenes)
        overlapMatrix.append(numOL)
    return overlapMatrix


def markJoanaActive(joanaOutput:list)->list:
    """
    helper function to mark JOANA-active terms

    Parameters
    ----------
    joanaOutput : list
        output of JOANA method.

    Returns
    -------
    active : list
        true/false list depending on whether gene set is marked active or not.

    """
    active = []
    for t in joanaOutput:
        active.append(True) if float(t[1]) >= .5 else active.append(False)      
    return active


def markPadogActive(res_padog:list)->list:
    """
    helper function to mark PADOG-active terms

    Parameters
    ----------
    res_padog : list
        output of PADOG method.

    Returns
    -------
    active : list
        true/false list depending on whether gene set is marked active or not.

    """
    active =[]
    for t in res_padog:
        active.append(True) if float(t[4]) < .05 else active.append(False)  
    return active


def markFisherActive(terms:list, asMatrix:list, joanaQValues:list)->list:
    """
    assess which terms are significantly enriched with the basic Fisher's exact test

    Parameters
    ----------
    terms : list
        description list of all terms.
    asMatrix : list
        assignment matrix of genes to gene sets.
    joanaQValues : list
        q-values for all genes.

    Returns
    -------
    active : list
        true/false list depending on whether a gene set is marked active or not.

    """
    active = []
    for i in range(0, len(terms)):
        oddsratio, pvalue = test_fisher(i, asMatrix, joanaQValues)
        active.append(True)  if pvalue <= .05 else active.append(False)

    return active
                 


def test_fisher(index:int, asMatrix:list, joanaQValues:list)->tuple:
    """
    function to identify gene sets that are enriched compared to the others gene sets

    Parameters
    ----------
    index : int
        index of the gene set to analyze.
    asMatrix : list
        assignment matrix.
    joanaQValues : list
        q-values.

    Returns
    -------
    oddsratio : float
        odds ratio of stats.fisher_exact.
    pvalue : float
        p-value of stats.fisher_exact.

    """
    qBiggerAnno = 0
    qSmallerAnno = 0
    qBiggerNotAnno = 0
    qSmallerNotAnno = 0
    
    for i in range(0, len(asMatrix)):

        qValue = joanaQValues[i]
        if index in asMatrix[i]:
            if qValue <= THRESHOLD:
                qSmallerAnno += 1
            else:
                qBiggerAnno += 1
        else:
            if qValue <= THRESHOLD:
                qSmallerNotAnno += 1
            else: 
                qBiggerNotAnno += 1
        
    oddsratio, pvalue = stats.fisher_exact([[qSmallerAnno, qSmallerNotAnno], [qBiggerAnno, qBiggerNotAnno]])
    return oddsratio, pvalue


def getMalRel(activeTerms1:list, activeTerms2:list, terms:list, malKegg:list)->list:
    """
    assess the two active Terms lists and if either term is active check relevance score and translate that score into color code

    Parameters
    ----------
    activeTerms1 : list
        list of active/inactive terms for method 1, e.g., JOANA.
    activeTerms2 : list
        list of active/inactive terms for method 2.
    terms : list
        description of all terms.
    malKegg : list
        malacards relevance scores.

    Returns
    -------
    colorMalacard : list
        list of rgb colors based on malacards relevance scores.

    """
    colorMalacard = []


    for i in range(0, len(activeTerms1)):
        if activeTerms1[i] or activeTerms2[i]:
            name = terms[i][:8]
            
            perc = -0.4
            color = get_rgb_color(perc)
            
            for j in range(0, len(malKegg)):
                if name == malKegg[j][0] or "\""+name+"\"" == malKegg[j][0]:
                    perc = float(malKegg[j][2])/100.0
                    color = get_rgb_color(perc)
            colorMalacard.append(color)
            
    return  colorMalacard


def singleTermAnalysis(term:str, actTerms1:list, actTerms2:list, allTerms:list, jqvalues:list, asMatrix:list)->list:
    """
    compares the overlap of one term of interest with all other active terms of the one or two methods.
    return colors code for active nodes to display number overlap of significantly enriched genes with colors

    Parameters
    ----------
    term : str
        term of interest.
    actTerms1 : list
        active terms of method 1.
    actTerms2 : list
        active terms of method 2.
    allTerms : list
        description of all terms.
    jqvalues : list
        joana q-values.
    asMatrix : list
        assignment matrix.

    Returns
    -------
    list
        color codes for the graphs of the active terms.

    """
    

    termDivi, _ = compareTwoTerms(term, term, allTerms, asMatrix, jqvalues)

    colors = []  
    for j in range(0, len(allTerms)):
        bothSignif = 0
        if actTerms1[j] or actTerms2[j]:
            term2ID = allTerms[j]
        
            bothSignif, bothNotSig = compareTwoTerms(term, term2ID, allTerms, asMatrix, jqvalues)
            col = get_rgb_color((bothSignif / termDivi)*4)
            colors.append(col)
    return colors


    
def ActiveOverlap(ovlMatrix:list, actTerms1:list, actTerms2:list, shape1:str, shape2:str, terms:list, malKegg:list, filename:str, singleTermName:str = None, jqvalues:list = None, asMatrix:list = None)->tuple:
    """
    starts the overlap analysis and gathers node and edge information needed to plot active-terms graphs

    Parameters
    ----------
    ovlMatrix : list
        overlap matrix of gene sets.
    actTerms1 : list
        active terms of method 1.
    actTerms2 : list
        active terms of method 2.
    shape1 : str
        node shape designated for method 1 terms.
    shape2 : str
        node shape designated for method 2 terms.
    terms : list
        description of all terms.
    malKegg : list
        malacards-kegg relevance score scores.
    filename : str
        name for plot file.
    singleTermName : TYPE, optional
        name for single term of interest to plot overlap of significant terms with it. The default is None:str.
    jqvalues : TYPE, optional
        joana q-values. The default is None:list.
    asMatrix : TYPE, optional
        assignment matrix. The default is None:list.

    Returns
    -------
    tuple
        returns node degree analysis results from node degree method.

    """
    names = []
    edges = []
    edgeSize = []
    nodeShape = []
    nodeSize = []
    jointShape = "circle" if actTerms1 != actTerms2 else shape1
    
        
    vertexCounterOuter = 0
    
    
    for i in range(0,len(ovlMatrix)):
        
        if actTerms1[i] and actTerms2[i]:
            nodeShape.append(jointShape)
        
        elif actTerms1[i]:
            nodeShape.append(shape1)
            
        elif actTerms2[i]:            
            nodeShape.append(shape2)
            
        if actTerms1[i] or actTerms2[i]:
            
            names.append(terms[i][:8])
            nodeSize.append(ovlMatrix[i][i]/5)
            
            vertexCounterInner = vertexCounterOuter+1
            for h in range(i, len(ovlMatrix[i])):
                if i!=h and (actTerms1[h] or actTerms2[h]) and ovlMatrix[i][h]>0:
                    edges.append((vertexCounterOuter, vertexCounterInner))
                    edgeSize.append(math.log(ovlMatrix[i][h]/10))
                    vertexCounterInner+=1

            vertexCounterOuter+=1
    
    if singleTermName is not None:
        malacardColorCode = singleTermAnalysis(singleTermName, actTerms1, actTerms2, terms, jqvalues, asMatrix)
    else:
        malacardColorCode = getMalRel(actTerms1, actTerms2, terms, malKegg)
    
    plotGraph(names, edges, edgeSize, nodeSize, nodeShape, malacardColorCode, filename)
    return node_Degree(edges)


def plotGraph(names:list, edges:list, edgeSize:list, nodeSize:list, nodeShape:list, nodeColor:list, filename:str):
    """
    function to plot the graphs the result from the analysis

    Parameters
    ----------
    names : list
        descriptor names of active terms.
    edges : list
        edges connected the terms/nodes.
    edgeSize : list
        thickness list of the edges displaying the number of overlapping genes.
    nodeSize : list
        node size displaying the term/gene set size, i.e., the number of genes annotated to it.
    nodeShape : list
        the node shape describes which method(s) has (have) activated the term.
    nodeColor : list
        node color to describe the malacards relevance score.
    filename : str
        filename to save plot under.

    Returns
    -------
    None.

    """
    
    tree = Graph(directed = False)
    tree.add_vertices(len(names))
    
    ed = []
    edsi = []
    for i in range(0, len(edges)):
        if edgeSize[i] >= 1 or True:
            ed.append(edges[i])
            edsi.append(edgeSize[i])
    
    tree.add_edges(ed)
    
    visual_style = {}
    visual_style["vertex_size"] = nodeSize
    visual_style["vertex_label"] = names
    visual_style["vertex_shape"] = nodeShape
    visual_style["edge_width"] = edsi
    visual_style["bbox"] = (1000, 1000)
    visual_style["margin"] = 50
    visual_style["vertex_label_color"] = ["red" for n in names]
    visual_style["vertex_color"] = nodeColor
    visual_style["vertex_label_dist"] = [1 for n in names]
    out = plot(tree,layout = tree.layout("kk"),  **visual_style)

    out.save(OUTPUT_PATH + filename + ".png")


def compareTwoTerms(term1:str, term2:str, allTerms:list, assignMatrix:list, jQValues:list, additional_info_text:str = None)->tuple:
    """
    overlap comparison between two terms of intereset, detecting how many genes overlap
    or are present in one of the two terms and are significantly enriched or not.

    Parameters
    ----------
    term1 : str
        term 1 to assess overlap with.
    term2 : str
        term 2 to assess overlap with.
    allTerms : list
        description of all terms.
    assignMatrix : list
        assignment matrix .
    jQValues : list
        q-values of JOANA.
    additional_info_text : TYPE, optional
        optional info text to be output in the console next to overlap analysis results. The default is None:str.

    Returns
    -------
    tuple
        DESCRIPTION.

    """
    ind1 = allTerms.index(term1)
    ind2 = allTerms.index(term2)
    bo = 0
    bono = 0
    t1 = 0
    t1no = 0
    t2 = 0
    t2no = 0
    
    for i in range(0, len(assignMatrix)):
        qval = jQValues[i]

        if ind1 in assignMatrix[i] and ind2 in assignMatrix[i]: #both term are annotated in gene
            if qval <= THRESHOLD:
                bo+=1
            else:
                bono +=1
        elif ind1 in assignMatrix[i]:
            if qval <= THRESHOLD:
                t1+=1
            else:
                t1no+=1
        elif ind2 in assignMatrix[i]:
            if qval <= THRESHOLD:
                t2+=1
            else:
                t2no+=1
    
    print("Genes significant in both,", bo, ", in both not, ", bono, "," , term1[:8] , ",", t1, t1no, "," , term2[:8], ",", t2, t2no, ",", additional_info_text)
    return bo, bono        
     

def analyzeDatasets():
    """
    Method to analyze a set of datasets at once.
    Skipping datasets that do not contain 338 gene sets or unsupported formatting
    Method produces graph for JOANA, PADOG and Fisher's exact test and their joint graphs for terms marked active by them
    Further anaylsis may be specified for certain datasets of interest.

    Returns
    -------
    None.

    """
    global FILEPATH_ASSIGNMENT_MATRIX, FILEPATG_KEGG, FILEPATH_MALACARD_KEGG, FILEPATH_QVALUES_GENEID
    global FILEPATH_RES_PADOG, FILEPATH_OUTPUT_JOANAPY, FILEPATH_QVALUES_JOANA, FILEPATH_TERMS, OUTPUT_PATH
    
    for i in dataset_list:
        print("Start", i)
        try:
            FILEPATH_ASSIGNMENT_MATRIX = basepath + i + "/assignment_matrix_filtered_min0_max100000.txt"
            FILEPATH_KEGG = basepath + "kegg.gmt"
            FILEPATH_MALACARD_KEGG = basepath + i + "/malacard_kegg.csv"
            FILEPATH_QVALUES_GENEID = basepath + i + "/" + i + "_qvalues.csv"
            FILEPATH_RES_PADOG = basepath + i + "/res_padog.csv"
            FILEPATH_OUTPUT_JOANAPY = basepath + i + "/output_joanapy_single_species.csv"
            FILEPATH_QVALUES_JOANA = basepath + i + "/qvalues_first_species.txt"
            FILEPATH_TERMS  = basepath + i + "/terms.txt"
            OUTPUT_PATH = basepath + i + "/"
            
            asMatrix, kegg, malKegg, qValuesGeneId, joanaOutput, joanaQValues, terms, res_padog = loadAllFiles()
            j = markJoanaActive(joanaOutput)
            p = markPadogActive(res_padog)
            s = markFisherActive(terms, asMatrix, joanaQValues)
        
            overlapMatrix = keggOverlap(kegg)
            
            
            print("# active with JOANA:", j.count(True), "/ # active with Padog:", p.count(True), "/ # active with Fisher:", s.count(True))
            print(len(j), len(p), len(s))
            
            jmean, jvar, jnum, jdegree = (0,0,0,0)
            pmean, pvar, pnum, pdegree = (0,0,0,0)
            if j.count(True) > 0:
                jmean, jvar, jnum, jdegree = ActiveOverlap(overlapMatrix, j, j, "square", "square", terms, malKegg, "JOANAoverlap")
            
            if s.count(True) > 0:
                ActiveOverlap(overlapMatrix, s, s, "triangle-up", "triangle-up", terms, malKegg, "Fisheroverlap")
        
            
            if j.count(True) > 0 and s.count(True) > 0:
                ActiveOverlap(overlapMatrix, j, s, "square", "triangle-up", terms, malKegg, "JOANAandFisher")
                
            if p.count(True) > 0 and j.count(True) > 0:
                ActiveOverlap(overlapMatrix, j, p, "square", "diamond", terms, malKegg, "JOANAandPADOG")
                if i == "GSE14924_CD4":
                    ActiveOverlap(overlapMatrix, j, p, "square", "diamond", terms, malKegg, "JOANAPADOGsingletermoverlapmap", "hsa05168_Herpes_simplex_virus_1_infection", joanaQValues, asMatrix)
                elif i == "GSE24739_G0":
                    ActiveOverlap(overlapMatrix, j, p, "square", "diamond", terms, malKegg, "JOANAPADOGsingletermoverlapmap", "hsa05206_MicroRNAs_in_cancer", joanaQValues, asMatrix)
            if p.count(True) > 0 and s.count(True) > 0:
                ActiveOverlap(overlapMatrix, p, s, "diamond", "triangle-up", terms, malKegg, "PADOGandFisher")
                
            if p.count(True) > 0:
                pmean, pvar, pnum, pdegree = ActiveOverlap(overlapMatrix, p, p, "diamond", "diamond", terms, malKegg, "PADOGoverlap")
            
            tval, pval = stats.ttest_ind(jdegree, pdegree, equal_var=True)
            
            print("scipy t_value", tval, "pval", pval)
            print("pmean", pmean, "pnum", pnum, "jmean", jmean, "jnum", jnum )
            
            if i == "GSE14924_CD4":
                print("CD4 term comparison")
                compareTwoTerms("hsa05168_Herpes_simplex_virus_1_infection", "hsa05131_Shigellosis", terms, asMatrix, joanaQValues, "Rather big Padog only term")
                compareTwoTerms("hsa05168_Herpes_simplex_virus_1_infection", "hsa05202_Transcriptional_misregulation_in_cancer", terms, asMatrix, joanaQValues, "Very relevant Padog only term")
                compareTwoTerms("hsa05168_Herpes_simplex_virus_1_infection", "hsa05322_Systemic_lupus_erythematosus", terms, asMatrix, joanaQValues, "medium size and relavant, active with both")
                compareTwoTerms("hsa05168_Herpes_simplex_virus_1_infection", "hsa05164_Influenza_A", terms, asMatrix, joanaQValues, "Joana only, fairly big and relevant")
                compareTwoTerms("hsa05168_Herpes_simplex_virus_1_infection", "hsa05203_Viral_carcinogenesis", terms, asMatrix, joanaQValues, "Another Padog only term")
            
            if i == "GSE24739_G0":
                print("G0 term comparison")
                compareTwoTerms("hsa05206_MicroRNAs_in_cancer", "hsa05220_Chronic_myeloid_leukemia", terms, asMatrix, joanaQValues, "Padog only, very relevant\n")
                compareTwoTerms("hsa05206_MicroRNAs_in_cancer", "hsa05214_Glioma", terms, asMatrix, joanaQValues, "Padog only, somewhat relevant\n")
                compareTwoTerms("hsa05206_MicroRNAs_in_cancer", "hsa05230_Central_carbon_metabolism_in_cancer", terms, asMatrix, joanaQValues, "Both active, somewhat relevant\n")    
                compareTwoTerms("hsa05206_MicroRNAs_in_cancer", "hsa05212_Pancreatic_cancer", terms, asMatrix, joanaQValues, "Both active, very relevant\n")
                compareTwoTerms("hsa05206_MicroRNAs_in_cancer", "hsa05166_Human_T-cell_leukemia_virus_1_infection", terms, asMatrix, joanaQValues, "Padog only, somewhat relevant, big set size\n")
                compareTwoTerms("hsa05206_MicroRNAs_in_cancer", "hsa05210_Colorectal_cancer", terms, asMatrix, joanaQValues, "JOANA only, very relevant\n")
                    
            
            sig = 0
            for q in joanaQValues:
                if q <= THRESHOLD:
                    sig += 1
            print("Assessment of all JOANA-qvalues in :", i)
            print("Below threshold:", sig, "total:", len(joanaQValues), "ratio:", sig/len(joanaQValues))
            
        except IndexError as e:
            print(e, "for gene set", i)
        except TypeError as e:
            print(e, "none return from node_Degree")            
        except UnboundLocalError as e:
            print(e)
        print("End", i)


def node_Degree(edges:list)->tuple:
    """
    extracts all node degree information through list of graph edges

    Parameters
    ----------
    edges : list
        list of edges in a single-method active-terms graph .

    Returns
    -------
    tuple
        tuple of mean, variance of node degrees and number of node and list of node degrees.

    """
    
    if len(edges) == 0:
        return
    d = {}

    for i in range(0, len(edges)):
        n1 = edges[i][0]
        n2 = edges[i][1]
        if n1 not in d:
            d[n1]=1
        else:
            d[n1]+=1
        if n2 not in d:
            d[n2]=1
        else:
            d[n2]+=1
    m = mean(list(d.values()))
    v = variance(list(d.values()))
    
    return m, v, edges[-1][1]+1, list(d.values())



            

def main():
    
    """
    extracted analysis for dataset GSE24739_G0
    

    asMatrix, kegg, malKegg, qValuesGeneId, joanaOutput, joanaQValues, terms, res_padog = loadAllFiles()
    
    j = markJoanaActive(joanaOutput)
    p = markPadogActive(res_padog)
    s = markFisherActive(terms, asMatrix, joanaQValues)
    
    overlapMatrix = keggOverlap(kegg)
    
    
    print("# active with JOANA:", j.count(True), "/ # active with Padog:", p.count(True), "/ # active with Fisher:", s.count(True))
    print(len(j), len(p), len(s))
    
    ActiveOverlap(overlapMatrix, j, j, "square", "square", terms, malKegg, "JOANAoverlap")
    ActiveOverlap(overlapMatrix, p, p, "diamond", "diamond", terms, malKegg, "PADOGoverlap")
    ActiveOverlap(overlapMatrix, s, s, "triangle-up", "triangle-up", terms, malKegg, "Fisheroverlap")

    ActiveOverlap(overlapMatrix, j, p, "square", "diamond", terms, malKegg, "JOANAandPADOG")
    ActiveOverlap(overlapMatrix, j, s, "square", "triangle-up", terms, malKegg, "JOANAandFisher")
    ActiveOverlap(overlapMatrix, p, s, "diamond", "triangle-up", terms, malKegg, "PADOGandFisher")
    
    """
    
    analyzeDatasets()
    
if __name__ == "__main__":
    main()
