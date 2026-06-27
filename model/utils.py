import math 

def mean_var(obj : list): 
    
    total = 0

    if len(obj) == 0:
        raise ValueError("invalid list length")

    for item in obj:
        total+=item
    
    mean = total/len(obj)

    variance = 0 
    for item in obj: 
        var_to_add = (item-mean)**2

        variance+= var_to_add
    
    variance = variance/len(obj)
    return mean,variance


def segement_score (obj: list, tau: int): 

    if tau>(len(obj)-1) or (tau == 0): 
        raise ValueError("invalid tau index")

    list_left = obj[:tau]
    list_right = obj[tau:]

    mean_1,var_1 = mean_var(list_left)
    mean_2,var_2 = mean_var(list_right)

    return  -(len(list_left)*var_1+len(list_right)*var_2)


def softmax(obj:dict): 
    scores = []
    
    scores = list(obj.values())

    max_num = max(scores)

    for i in range(0,len(scores)): 
        scores[i] = scores[i]-max_num
    
    for i in range (0,len(scores)): 
        scores[i] = math.exp(scores[i])

    total = 0 
    for i in range(0,len(scores)): 
        total+=scores[i]

    for i in range(0,len(obj)): 
        obj[i+2] = scores[i]/total
    
    return obj 

if __name__ == "__main__": 
    print(mean_var([2,4,6]))