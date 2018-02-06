#%matplotlib inline
import statemap
import matplotlib.pyplot as plt
import numpy as np

file = "./president_general_polls_2016.csv"


poll_data = statemap.election.poll(file)

#file = "./test.csv"
#e = election(file)

def bar_plot(quantity, trump, clinton):
    '''taken from stack overflow thread about horizontal bar charts and modified for this class
    Inputs: quantity used differentiate between electoral and popular vote; trump and clinton represent data for each candidate
    '''
    #quantity = ('%')
    segments = 2
    if quantity == ("E"):
        strs = [' votes', 'Electoral Votes']
    elif quantity == ("%"):
        strs = ['% ', 'Popular Votes']
        
    # generate some multi-dimensional data & arbitrary labels
    data =(np.array(([trump, clinton])))
    total = (np.array(([trump, clinton])))
    y_pos = np.arange(len(quantity))
    fig = plt.figure(figsize=(12,1))
    ax = fig.add_subplot(111)
    candidates = ["Trump","Clinton"]
    colors ='rb'
    patch_handles = []
    left = np.zeros(1) # left alignment of data starts at zero
    for i, d in enumerate(data):
        patch_handles.append(ax.barh(y_pos, d, 
          color=colors[i%len(colors)], align='center', 
          left=left))
        # accumulate the left-hand offsets
        left += d

    # go through all of the bar segments and annotate
    i = 0
    for j in range(len(patch_handles)):
        for i, patch in enumerate(patch_handles[j].get_children()):
            bl = patch.get_xy()
            x = 0.5*patch.get_width() + bl[0]
            y = 0.5*patch.get_height() + bl[1]

            ax.text(x,y, "{}: {}{}".format(candidates[j],total[j],strs[0]), ha='center')


    ax.set_yticks(y_pos)
    ax.set_yticklabels(quantity)
    ax.set_xlabel(strs[1])

    plt.show()



def popular_vote(p):

    trump_pop_votes = 0
    clinton_pop_votes = 0
    
    for s in p:

        clinton = p[s][0]
        trump = p[s][1]

        electoral_votes = statemap.election.electorate()[s]
        clinton_result = round(np.random.normal(loc=clinton),2)
        trump_result = round(np.random.normal(loc=trump),2)
        if clinton != "nan" and trump != "nan":
            clinton_pop_votes += clinton_result*float((electoral_votes/538))
            trump_pop_votes += trump_result*float((electoral_votes/538))

    clinton_pop_votes = round(clinton_pop_votes,2)
    trump_pop_votes = round(trump_pop_votes,2)
    
    bar_plot(("%"), trump_pop_votes, clinton_pop_votes)

    

def electoral_vote(p):
    trump_electoral_votes = 0
    clinton_electoral_votes = 0
    
    for s in p:

        clinton = p[s][0]
        trump = p[s][1]

        electoral_votes = statemap.election.electorate()[s]
        clinton_result = clinton
        trump_result = trump
        #clinton_result = round(np.random.normal(scale=sd,loc=clinton),2)
        #trump_result = round(np.random.normal(scale=sd,loc=trump),2)
        
        if clinton_result > trump_result:
            clinton_electoral_votes += electoral_votes
        elif trump_result > clinton_result:
            trump_electoral_votes += electoral_votes
        else:
            pass
        
    bar_plot(("E"), trump_electoral_votes, clinton_electoral_votes)
    
#The code below is left unused but served as a template for how to calculate votes or show the visualization
#Keeping it in case my code above messes up

def calc_votes(p):
    trump_electoral_votes = 0
    clinton_electoral_votes = 0
    trump_pop_votes = 0
    clinton_pop_votes = 0
    
    for s in p:

        clinton = p[s][0]
        trump = p[s][1]

        electoral_votes = statemap.election.electorate()[s]
        clinton_result = round(np.random.normal(scale=sd,loc=clinton),2)
        trump_result = round(np.random.normal(scale=sd,loc=trump),2)
        clinton_pop_votes += clinton_result*float((electoral_votes/538))
        trump_pop_votes += trump_result*float((electoral_votes/538))

        if clinton_result > trump_result:
            clinton_electoral_votes += electoral_votes
        elif trump_result > clinton_result:
            trump_electoral_votes += electoral_votes
        else:
            pass
            #print("No result reported")
    clinton_pop_votes = round(clinton_pop_votes,2)
    trump_pop_votes = round(trump_pop_votes,2)
    result_string = '''Election Results:
    Popular vote - Trump: {}% Clinton: {}%
    Electoral vote - Trump: {} Clinton: {}'''.format(trump_pop_votes,clinton_pop_votes,trump_electoral_votes,clinton_electoral_votes)
    
    return (trump_electoral_votes,clinton_electoral_votes,trump_pop_votes,clinton_pop_votes)



def default_plot():
    
    p = statemap.election.poll()
    for s in p:

        clinton = p[s][0]
        trump = p[s][1]

        electoral_votes = statemap.election.electorate()[s]
        clinton_result = round(np.random.normal(scale=sd,loc=clinton),2)
        trump_result = round(np.random.normal(scale=sd,loc=trump),2)
        clinton_pop_votes += clinton_result*float((electoral_votes/538))
        trump_pop_votes += trump_result*float((electoral_votes/538))

        if clinton_result > trump_result:
            clinton_electoral_votes += electoral_votes
        elif trump_result > clinton_result:
            trump_electoral_votes += electoral_votes
        else:
            pass
            #print("No result reported")
    clinton_pop_votes = round(clinton_pop_votes,2)
    trump_pop_votes = round(trump_pop_votes,2)
    result_string = '''Election Results:
    Popular vote - Trump: {}% Clinton: {}%
    Electoral vote - Trump: {} Clinton: {}'''.format(trump_pop_votes,clinton_pop_votes,trump_electoral_votes,clinton_electoral_votes)

    percent = ('%')
    segments = 2

    # generate some multi-dimensional data & arbitrary labels
    data =(np.array(([trump, clinton])))
    percentage = (np.array(([trump, clinton])))
    y_pos = np.arange(len(percent))

    fig = plt.figure(figsize=(8,1))
    ax = fig.add_subplot(111)
    candidates = ["Trump","Clinton"]
    colors ='rb'
    patch_handles = []
    left = np.zeros(len(percent)) # left alignment of data starts at zero
    for i, d in enumerate(data):
        patch_handles.append(ax.barh(y_pos, d, 
                                     color=colors[i%len(colors)], align='center', 
                                     left=left))
        # accumulate the left-hand offsets
        left += d

        # go through all of the bar segments and annotate
        i = 0
        for j in range(len(patch_handles)):
            for i, patch in enumerate(patch_handles[j].get_children()):
                bl = patch.get_xy()
                x = 0.5*patch.get_width() + bl[0]
                y = 0.5*patch.get_height() + bl[1]

                ax.text(x,y, "%d%%" % (percentage[j]) + ' ' + candidates[j], ha='center')


                ax.set_yticks(y_pos)
                ax.set_yticklabels(percent)
                ax.set_xlabel('Popular Vote Percentage')
                ax.set_xticks([])

    plt.show()

#vote_tup = calc_votes(poll_data)
#t_ev = vote_tup[0]
#c_ev = vote_tup[1]
#t_pv = vote_tup[2]
#c_pv = vote_tup[3]