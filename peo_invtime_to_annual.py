#function to get recaccurance rate from poe deagg output

def poe_invtime_to_annual(poes, investigation_time):
    from numpy import array, log
    # now get annualised curves
    P0 = 1 - array(poes)
    n = -1*log(P0)
    annual_poes = n / investigation_time
    return annual_poes

print(1/poe_invtime_to_annual(0.00110349, 50))

