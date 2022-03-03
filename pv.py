import math
from scipy.special import comb
import sys
import bisect

import main as mains


def p_upper_exact(n,N,M,Q,T):
    """
    according to the parameters, we choose the way to calculate p-value (the upper bound way: cal_pvalue_upper;  the exact way: cal_pvalue_simple)

    :param n: the number of internal nodes of ns
    :param N: the total nodes of the whole graph
    :param M: the total edge-number of the whole graph
    :param Q: the minimal number of internal edges
    :param T: the maximal number of external edges
    :return:
    """

    Nn = int(n * (N - n))
    n2 = int(n*(n-1)/2)
    Nn2 = int((N-n)*(N-n-1)/2)
    # penalty = mains.get_log_factorial(N,n)
    low2 = 0
    upper2 = T
    low1 = min(max(M-Nn2,min(n2,Q)),n2)
    upper1 = min(M-upper2,n2)
    if M-low1-upper2 > Nn2:
        return float('-inf')
    '''this is the denominator of the y_{Q+1,l2}^{in} and g_{n(N-n)-T+1}^{out} in the main paper when using the upper bound p-value'''
    qout1 = int(Nn - upper2 + 1)
    qout2 = int(M - low1 - upper2 + 1)
    qin1 = int(Nn2 - M + low1 + +low2 + 1)

    '''if the the denominator is 0, we use the exact pvalue formulation'''
    if qout1 == 0 or qout2 == 0 or qin1 == 0:
        logp = cal_pvalue_simple(n, N, M,  Q, T)
    else:
        '''if meet the condition of y_{Q+1,l2}^{in} and g_{n(N-n)-T+1}^{out}, we use the upper bound pvalue '''
        qout = (upper2 * (Nn2 - M + low1 + upper2) * 1.0) / (qout1 * qout2 * 1.0)
        qin = ((n2 - low1) * (M - low1 - low2) * 1.0 / ((low1 + 1) * qin1 * 1.0))
        if qin < 1 and qout < 1 and qin >=0 and qout >= 0:
            logp = cal_pvalue_upper(n, N, M, low1, low2, upper1, upper2, qout, qin)
        else:
            logp = cal_pvalue_simple(n,N,M,Q,T)
    return logp

def cal_pvalue_upper(n,N,M,low1,low2,upper1,upper2,qout,qin):
    """
    according to the formulate [12] of the paper MSSC to calculate the p-value

    :param n:
    :param N:
    :param M:
    :param low1:
    :param low2:
    :param upper1:
    :param upper2:
    :param qout:
    :param qin:
    :return:
    """
    if n == 1:
        return 0
    Nn = int(n*(N-n));n2 = int(n*(n-1)/2);N2= int(N*(N-1)/2);Nn2= int((N-n)*(N-n-1)/2)
    # temp0 = mi.get_log_factorial(int(N2)) - mi.get_log_factorial(int(N2 - M)) - mi.get_log_factorial(int(M))
    # temp1 = mi.get_log_factorial(int(Nn2)) - mi.get_log_factorial(int(M-low1-upper2)) - mi.get_log_factorial(int(Nn2 - M + low1+upper2))
    # temp5 = mi.get_log_factorial(int(n2)) - mi.get_log_factorial(int(low1)) - mi.get_log_factorial(int(n2-low1))
    # temp2 = mi.get_log_factorial(int(Nn)) - mi.get_log_factorial(int(upper2)) - mi.get_log_factorial(int(Nn - upper2))
    temp0 = mains.get_log_factorial(N2, M)  # factorial_num[int(N2)] - factorial_num[int(N2 - M)] - factorial_num[int(M)]
    temp1 = mains.get_log_factorial(Nn2,M - low1 - upper2)  # factorial_num[int(Nn2)] - factorial_num[int(M-low1-upper2)] - factorial_num[int(Nn2 - M + low1+upper2)]
    temp5 = mains.get_log_factorial(n2,low1)  # factorial_num[int(n2)] - factorial_num[int(low1)] - factorial_num[int(n2-low1)]
    temp2 = mains.get_log_factorial(Nn,upper2)  # factorial_num[int(Nn)] - factorial_num[int(upper2)] - factorial_num[int(Nn - upper2)]

    temp3 = math.log((1 - math.pow(qout, upper2 - low2 + 1)) / (1 - qout))
    temp4 = math.log((1-math.pow(qin,(upper1-low1+1)))/(1-qin))
    logp = temp1+temp2+temp3+temp4+temp5-temp0
    return logp


def cal_pvalue_simple(n,N,M,Q,T):
    """
    according to the formulate [6] of the paper MSSC to calculate the p-value, due to the original formulate is a1+a2+...+an, we use the log value of p, and thus log(a1+a2+...+an),
    to calculate the log value we use the method of   cal_log_add(a,b)

    :param n:
    :param N:
    :param M:
    :param Q:
    :param T:
    :return:
    """
    calculate_time  = 0
    if n == 1:
        return 0
    Nn = int(n*(N-n));n2 = int(n*(n-1)/2);N2= int(N*(N-1)/2);Nn2= int((N-n)*(N-n-1)/2)
    low1 = int(min(n2,Q))
    upper1 = int(min(M, n2))
    stime = 0
    logp_exact = []
    break2 = False
    upper2 = int(min(Nn,min(M-low1,T)))
    if M-low1-upper2 > Nn2:
        return  float('-inf')
    # upper2 = int(min(Nn, min(M - low1-Nn2, T)))
    # upper2 = int(alpha_for_edge*(min(u2_1,n_out))+(1-alpha_for_edge)*(min(gupper,u2_1)))
    # temp_4 = mi.get_log_factorial(int(N2)) - mi.get_log_factorial(int(N2 - M)) - mi.get_log_factorial(int(M))
    # temp_1 = mi.get_log_factorial(int(n2)) - mi.get_log_factorial(int(low1)) - mi.get_log_factorial(int(n2 - low1))
    # temp_2 = mi.get_log_factorial(int(Nn)) - mi.get_log_factorial(int(upper2)) - mi.get_log_factorial(int(Nn - upper2))
    # temp_3 = mi.get_log_factorial(int(Nn2)) - mi.get_log_factorial(int(M - low1 - upper2)) - mi.get_log_factorial(int(
    #     Nn2 - M + low1 + upper2))
    temp_4 = mains.get_log_factorial(N2, M)  # factorial_num[int(N2)] - factorial_num[int(N2 - M)] - factorial_num[int(M)]
    temp_1 = mains.get_log_factorial(n2, low1)  # factorial_num[int(n2)] - factorial_num[int(low1)] - factorial_num[int(n2 - low1)]
    temp_2 = mains.get_log_factorial(Nn,upper2)  # factorial_num[int(Nn)] - factorial_num[int(upper2)] - factorial_num[int(Nn - upper2)]
    temp_3 = mains.get_log_factorial(Nn2, (M - low1 - upper2))

    temp_exact = temp_1 + temp_2 + temp_3 - temp_4
    logp_exact.append(temp_exact)
    '''i is the range of the internal edge number'''
    for i in range(low1, upper1+1):
        low2 = int(min(Nn,max(0,M-i-Nn2)))
        upper2 = int(min(min(M-i,T),Nn))
        # upper2 = int(min(min(M - i-Nn2, T), Nn))
        #temp_1 = mi.get_log_factorial(int(n2)) - mi.get_log_factorial(int(i)) - mi.get_log_factorial(int(n2 - i))
        temp_1 = mains.get_log_factorial(n2, i)
        if i == low1:
            j = upper2 - 1
        else:
            j = upper2
        while j >= low2:
            calculate_time += 1
            stime += 1
            # temp_2 = mi.get_log_factorial(int(Nn)) - mi.get_log_factorial(int(j)) - mi.get_log_factorial(int(Nn - j))
            # temp_3 = mi.get_log_factorial(int(Nn2)) - mi.get_log_factorial(int(M - i - j)) - mi.get_log_factorial(int(
            #     Nn2 - M + i + j))
            temp_2 = mains.get_log_factorial(Nn,j)
            temp_3 = mains.get_log_factorial(Nn2, M - i - j)
            temp_exact = temp_1+temp_2+temp_3-temp_4
            temp_exact_p = cal_log_add(temp_exact,logp_exact[stime-1])
            logp_exact.append(temp_exact_p)
            '''if the difference between the new log value of ai and the log(a1+..+ai-1) is bigger than 200, we breank the process'''
            if math.fabs(temp_exact-logp_exact[stime-1]) > 200:
                break2 = True
                break
            j=j-1
        if break2 == True:
            break

    return logp_exact[len(logp_exact)-1]

def cal_Q_g(W_inc,W_de,W_in,W_out):
    """
    according to the total internal edge-weight, we binary search the minimal number of edges whose total edge-weight is not smaller than W_in, and binary search the maximal number of edges whose total edge-weight is not bigger than W_out

    :param W_inc: the sum list of the increased weight order . the increase order of the normalized edge weight is {w1, w2, ..., wM} where, w1 <= w2 <= ... <= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...+wM}
    :param W_de: similar to W_inc, the sum list of the decreased weight order . the decrease order of the normalized edge weight is {w1, w2, ..., wM} where, w1 >= w2 >= ... >= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...wM}
    :param W_in: the total internal edge-weight
    :param W_out: the total external edge-weight
    :return: the minimal internal edge-number restriction due to the internal edge-weight, and the maximal external edge-number restriction due to the external edge-weight
    """
    length = len(W_inc)
    tempQ = bisect.bisect_right(W_inc, W_in, 0, length-1)
    if W_in < W_inc[0]:
        Qlow = 1
    elif W_inc[tempQ - 1] == W_in:
        Qlow = tempQ
    else:
        Qlow = tempQ + 1

    tempg = bisect.bisect_left(W_de, W_out, 0, length-1)
    if W_out < W_de[0]:
        gupper = 0
    elif W_inc[tempg] == W_in:
        gupper = tempg + 1
    else:
        gupper = tempg
    return Qlow, gupper

def cal_log_add(a,b):
    """
    calculate log(a+b) with loga and logb

    :param a:
    :param b:
    :return:
    """
    if a <= b:
        return b + math.log(1+math.exp(a-b))
    else:
        return a + math.log(1+math.exp(b-a))


''''
def binary_search_Q(low_index,upper_index,W_list,target,flag):
    #print low_index,upper_index,target
    li = int(low_index);ui=int(upper_index)
    if li == 0:
        if target <= W_list[0]:
            return 1
    if flag == 0:  # target > winbefore
        if W_list[li] >= target:
            return li + 1
        elif W_list[li+1] >= target:
            return li +2
    else:#target < winbefore <= Wlist[ui]
        if W_list[ui - 1] < target:
            return ui + 1
    while li < ui:
        mid = int(math.ceil((float(li) + ui)/2))
        value = W_list[mid]
        value_be = W_list[mid - 1]
        if value >= target and value_be < target:
            return mid + 1
        elif value < target:
            li = mid
        elif value_be > target:
            ui = mid-1
        else:
            return mid
    if li == ui:
            return li+1
def binary_search_g(low_index, upper_index, W_list, target,flag):
    # print "g:::",low_index, upper_index, target
    li = int(low_index)
    ui = int(upper_index)
    if li == 0:
        if target == W_list[0]:
            return 1
        elif target < W_list[0]:
            return 0
    if flag == 0:# target > woutbefore > Wlist[li]
        if W_list[li + 1] > target:
            return li+1
        elif W_list[li+1]==target:
            return li+2
    else:#target <woutbefore
        if W_list[ui] <= target:
            return ui + 1
        elif W_list[ui-1] <= target:
            return ui
    while li < ui:
        mid = int(math.ceil((float(li) + ui)/2))
        value = W_list[mid]
        value_be = W_list[mid + 1]
        if value <= target and value_be > target:
            return mid + 1
        elif value_be < target:
            li = mid
        elif value > target:
            ui = mid-1
        else:
            return mid + 2
    if li == ui:
            return li + 1
def cal_Q_g_binary(W_inc,W_de,W_in,W_out,QgWinWoutbefore):
    Qlbefore = QgWinWoutbefore[0]
    gubefore = QgWinWoutbefore[1]
    Winbefore = QgWinWoutbefore[2]
    Woutbefore = QgWinWoutbefore[3]
    length = len(W_inc)
    gupper = 0;Qlow = 1
    if W_in == Winbefore:
        Qlow = Qlbefore
    else:
        if W_in != 0:
            if W_in > Winbefore:
                # Qlow =
                Qlow = binary_search_Q(Qlbefore-1,length-1, W_inc, W_in,0)
            else:
                #flag = 1
                Qlow = binary_search_Q(0, Qlbefore-1, W_inc, W_in,1)

    if W_out == Woutbefore:
        gupper = gubefore
    else:
        if W_out != 0:
            if W_out > Woutbefore:
                gupper = binary_search_g(gubefore, length-1, W_de, W_out,0)
            else:
                gupper = binary_search_g(0, gubefore, W_de, W_out,1)
    QgWinWoutbefore[0] = Qlow
    QgWinWoutbefore[1] = gupper
    QgWinWoutbefore[2] = W_in
    QgWinWoutbefore[3] = W_out
    return Qlow, gupper,QgWinWoutbefore
'''