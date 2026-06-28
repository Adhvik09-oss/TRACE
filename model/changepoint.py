from utils import softmax,segement_score,mean_var
import math 
import statistics 

def score_splits(obj:list): 
    list_length = len(obj)

    scores = {}

    for tau in range(2,list_length-1): 


        scores[f"tau {tau}"] = segement_score(obj,tau)

    return scores

def fit_patient(ca:list,kras:list): 

    #declare floor values for ca and ct 

    ca_floor = math.log(37)
    ct_floor = math.log(5.0)

    #filter out for matching/  non nan values in both ca_alus and ct_values    
    valid = [i for i in range(min(len(ca), len(kras))) 
         if not (math.isnan(ca[i]) or math.isnan(kras[i]))]
    ca_values = [ca[i] for i in valid]
    ct_values = [kras[i] for i in valid]

    #now we have too look if the scores are valid for the patient 
    # if 80% or more of the values are floor values or the std is 

    non_informative_ca =0 
    non_informative_ct = 0

    # these are the calculations to see wheter or not all of our CA19-9/KRAS information is valid or not to use
    ca_informative = True
    ct_informative = True 

    for i in range(len(ca_values)): 
        if ca_values[i] <= ca_floor:
            non_informative_ca+=1
        if ct_values[i] <= ct_floor: 
            non_informative_ct+=1

    floor_frac_ca = non_informative_ca/len(ca_values)
    floor_frac_ct = non_informative_ct/len(ct_values)

    std_ca = statistics.stdev(ca_values)
    std_ct = statistics.stdev(ct_values)

    ca_informative = (floor_frac_ca<0.8) and (std_ca>0.3)
    ct_informative = (floor_frac_ct<0.8) and (std_ct>0.3)

    # only do 
    if ca_informative:

        ca19_scores = score_splits(ca_values)
        ca19_probs = softmax(ca19_scores)

    if ct_informative:    
        kras_scores = score_splits(ct_values)
        kras_probs = softmax(kras_scores)

    

    if ca_informative and ct_informative:
        mode = "dual"
    elif ca_informative:
        mode ="CA only"
        final_probs = ca19_probs
    elif ct_informative: 
        mode = "KRAS only"
        final_probs = kras_probs
    else: 
        return None 


    if mode == "dual":

        final_probs={}

        for key,val in ca19_probs.items(): 
            final_probs[key] = (ca19_probs[key]*kras_probs[key])
    
        total = sum(final_probs.values())
        final_probs = {k: v/total for k, v in final_probs.items()}

        best_tau = max(final_probs, key=final_probs.get)
        certainty = final_probs[best_tau]

        mu0_ca = mean_var(ca_values[:best_tau])[0]
        mu1_ca = mean_var(ca_values[best_tau:])[0]
        mu0_ct = mean_var(ct_values[:best_tau])[0]
        mu1_ct = mean_var(ct_values[best_tau:])[0]

        return {
    'tau': best_tau,
    'certainty': certainty,
    'probs': final_probs,
    'mu0_ca': mu0_ca,
    'mu1_ca': mu1_ca,
    'mu0_ct': mu0_ct,
    'mu1_ct': mu1_ct,
    'mode': mode
}
    elif mode == "CA only":

        best_tau = max(final_probs, key=final_probs.get)
        certainty = final_probs[best_tau]
        
        mu0_ca = mean_var(ca_values[:best_tau])[0]
        mu1_ca = mean_var(ca_values[best_tau:])[0]
        mu0_ct = None
        mu1_ct = None

        return {
    'tau': best_tau,
    'certainty': certainty,
    'probs': final_probs,
    'mu0_ca': mu0_ca,
    'mu1_ca': mu1_ca,
    'mu0_ct': mu0_ct,
    'mu1_ct': mu1_ct,
    'mode': mode
} 
    elif mode == "KRAS only":
        best_tau = max(final_probs, key=final_probs.get)
        certainty = final_probs[best_tau]

        mu0_ca = None
        mu1_ca = None
        mu0_ct = mean_var(ct_values[:best_tau])[0]
        mu1_ct = mean_var(ct_values[best_tau:])[0]


        return {
    'tau': best_tau,
    'certainty': certainty,
    'probs': final_probs,
    'mu0_ca': mu0_ca,
    'mu1_ca': mu1_ca,
    'mu0_ct': mu0_ct,
    'mu1_ct': mu1_ct,
    'mode': mode
}
        


 
                