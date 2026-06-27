from utils import softmax,segement_score


def score_splits(obj:list): 
    list_length = len(obj)

    scores = {}

    for tau in range(2,list_length-1): 


        scores[f"tau {tau}"] = segement_score(obj,tau)

    return scores

def fit_patient(obj:list): 
    
    scores =  score_splits(obj)
    scores = softmax(scores)

    taus = list(scores.values())
    max_tau = max(taus)

    for key,val in scores.items(): 
        if val == max_tau: 
            return key,val


        