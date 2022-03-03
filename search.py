import pv
import random
import global_var_model as gl
import math
def search_community(nb_dict,ns,W_ns_in, number_ns_in, W_ns_out, number_ns_out,edge_dict,link_table,N,M,W_inc,W_de,S):
    """
    firstly, we calculate the p-value of the candidate community ns, and then try to add every neighbor to ns ; try to remove every node from ns, and then choose a operation who has the most smaller p-value and update the ns

    :param nb_dict: the neighbors dict (whose key is the node index and value is the number of neightbors in ns)
    :param ns: the candidate community node set
    :param W_ns_in: the total internal edge-weight
    :param number_ns_in: the total number of internal edges
    :param W_ns_out: the total external edge-weight
    :param number_ns_out: the total number of external edges
    :param edge_dict: e.g. {node1: {node2: {weight: value1}, {node3: {weight: value2 }}}}, node1 and node2 is connected with edge-weight value1, node1 and node3 is connnected with edge-weight value2
    :param link_table: e.g.{node1: [node2,node3]}, the key is the node index and the value is the linked nodes list of the key node
    :param N: the total number of nodes
    :param M: the total number of edges
    :param W_inc: the sum list of the increased weight order . the increase order of the normalized edge weight is {w1, w2, ..., wM} where, w1 <= w2 <= ... <= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...+wM}
    :param W_de: similar to W_inc, the sum list of the decreased weight order . the decrease order of the normalized edge weight is {w1, w2, ..., wM} where, w1 >= w2 >= ... >= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...+wM}
    :param S: the node strength value dict, key is node index and value is node strength
    :param index_name: index and name reflact dict, key is node index and value is the original name
    :return: the result community and its p-value
    """

    wi = W_ns_in
    wo = W_ns_out
    ni = number_ns_in
    no = number_ns_out
    n= len(ns)
    Qlow, gupper = pv.cal_Q_g(W_inc, W_de,W_ns_in, W_ns_out)
    '''according to the parameter of weight influence and edge-number influence adjust the exact edge-number'''
    T = math.floor(gl.alpha_for_edge*number_ns_out+(1-gl.alpha_for_edge)*(gupper))
    Q = math.ceil(gl.alpha_for_edge*number_ns_in+(1-gl.alpha_for_edge)*(min(Qlow,n*(n-1)/2)))
    pvalue = pv.p_upper_exact(n,N,M,Q,T)
    while True:
        if pvalue == float('-inf'):
            return ns,pvalue
        na_ns = []
        for a in ns:
            na_ns.append(gl.index_name[a])
        add_node, min_add_p, awi, awo, ani, ano = add_function(nb_dict, pvalue, ns, N, M, W_inc, W_de, S, edge_dict,wi, wo, ni, no)
        remove_node,min_remove_p,rwi,rwo,rni,rno = remove_function(pvalue,ns,N,M,W_inc,W_de,S,edge_dict,wi, wo, ni, no,min_add_p)
        if gl.item >= gl.max_iter:
            return ns,pvalue
        if remove_node== -1 and add_node == -1:
            return ns,pvalue
        elif remove_node== -1:
            add_search_next(add_node,link_table,nb_dict,ns)
            pvalue = min_add_p
            wi = awi
            wo = awo
            ni = ani
            no = ano
        elif add_node == -1:
            remove_search_next(remove_node,link_table,nb_dict,ns)
            pvalue = min_remove_p
            wi = rwi
            wo = rwo
            ni = rni
            no = rno
        else:
            if min_remove_p < min_add_p:
                remove_search_next(remove_node, link_table, nb_dict, ns)
                pvalue = min_remove_p
                wi = rwi
                wo = rwo
                ni = rni
                no = rno
            elif min_remove_p > min_add_p:
                add_search_next(add_node, edge_dict, nb_dict, ns)
                pvalue = min_add_p
                wi = awi
                wo = awo
                ni = ani
                no = ano
            else:
                seed = random.randint(0,1)
                if seed == 0:
                    add_search_next(add_node,link_table,nb_dict,ns)
                    pvalue = min_add_p
                    wi = awi
                    wo = awo
                    ni = ani
                    no = ano
                else:
                    remove_search_next(remove_node,link_table,nb_dict,ns)
                    pvalue = min_remove_p
                    wi = rwi
                    wo = rwo
                    ni = rni
                    no = rno


def add_function(neib_dict, pvalue, ns, N, M, W_inc, W_de, S, edge_dict, W_ns_in, W_ns_out, number_ns_in,
                 number_ns_out):
    """
    for each neighbor node of ns, we try to find one node that added to ns and make the p-value of the new ns is much lower or not bigger than a threshold

    :param neib_dict: the neighbors dict (whose key is the node index and value is the number of neightbors in ns)
    :param pvalue_orig: the p-value of the original ns
    :param ns: the candidate community node set
    :param N: total node number
    :param M: total edge number
    :param W_inc:  the sum list of the increased weight order . the increase order of the normalized edge weight is {w1, w2, ..., wM} where, w1 <= w2 <= ... <= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...+wM}
    :param W_de: similar to W_inc, the sum list of the decreased weight order . the decrease order of the normalized edge weight is {w1, w2, ..., wM} where, w1 >= w2 >= ... >= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...+wM}
    :param S: the node strength value dict, key is node index and value is node strength
    :param edge_dict:  e.g. {node1: {node2: {weight: value1}, {node3: {weight: value2 }}}}, node1 and node2 is connected with edge-weight value1, node1 and node3 is connnected with edge-weight value2
    :param W_ns_in: the total internal edge-weight
    :param W_ns_out: the total external edge-weight
    :param number_ns_in: the total internal edge number
    :param number_ns_out: the total external edge number
    :return: the node can be add to ns; the p-value of new ns with the node added to ns; wi(wo): is the internal(external) edge-weight of the new ns;  ni(no) is the internal(external) edge number of new ns
    """

    wi = W_ns_in
    wo = W_ns_out
    ni = number_ns_in
    no = number_ns_out
    add_node = -1
    lst = []
    n = len(ns) + 1
    oldp = pvalue

    if len(neib_dict) == 0:
        return add_node, pvalue, wi, wo, ni, no
    # --------------------getting the informaitons of all the candidate add nodes in the neighbors of ns
    node_information = {}
    nin_nout = {}
    nin_win = {}
    nin_wout = {}

    for i in neib_dict:
        '''according to the original edge number, edge-weight values  and the ones connected with node i to adjuste the according value after remove i from ns'''
        Wns_node, WN_ns_node, number_ns_node, number_N_ns_node = cal_node_Min_Mout(ns, i, edge_dict, S)
        iin = W_ns_in + Wns_node
        iout = W_ns_out + WN_ns_node - Wns_node
        number_iin = number_ns_in + number_ns_node
        number_iout = number_ns_out + number_N_ns_node - number_ns_node
        '''store the win, wout, niin, nout information of the new ns (after add node i to ns) to node_information dict'''
        node_information[i] = (iin, iout, number_iin, number_iout)
        '''for two new candidate community a1(ns + node1), a2(ns + node2), if they have the same internal edge-number, we just need to calculate the one who has smaller external edge-number'''
        if number_iin in nin_nout:
            t1 = nin_nout[number_iin]
            if number_iout < t1[0]:
                nin_nout[number_iin] = [number_iout, i]
        else:
            nin_nout[number_iin] = [number_iout, i]
        '''for two new candidate community a1(ns + node1), a2(ns + node2), if they have the same internal edge-number, we just need to calculate the one who has bigger internal edge-weight'''
        if number_iin in nin_win:
            if iin > nin_win[number_iin][0]:
                nin_win[number_iin] = [iin, i]
        else:
            nin_win[number_iin] = [iin, i]
        '''for two new candidate community a1(ns + node1), a2(ns + node2), if they have the same internal edge-number, we just need to calculate the one who has smaller external edge-weight'''
        if number_iin in nin_wout:
            if iout < nin_wout[number_iin][0]:
                nin_wout[number_iin] = [iout, i]
        else:
            nin_wout[number_iin] = [iout, i]

    nodes_tobe_calculated = set()
    for nin in nin_nout:
        nodes_tobe_calculated.add(nin_nout[nin][1])
        nodes_tobe_calculated.add(nin_win[nin][1])
        nodes_tobe_calculated.add(nin_wout[nin][1])
    number_mout_minmum = float('inf')
    '''calculating the candidate nodes: if they have the same Q, we calculate the node who has the smallest T;
    if they have different Q,we calculate the node whose T is smaller than those Ts that have been touched
    until we meet the smallest T'''
    for node in nodes_tobe_calculated:
        temp = node_information[node]
        Qlow, gupper = pv.cal_Q_g(W_inc, W_de, temp[0], temp[1])
        Q = math.ceil(temp[2] * gl.alpha_for_edge + (min(Qlow, n * (n - 1) / 2)) * (1 - gl.alpha_for_edge))
        T = math.floor(temp[3] * gl.alpha_for_edge + gupper * (1 - gl.alpha_for_edge))
        number_mout_minmum = min(number_mout_minmum, T)
        # lst(minimual internal edge-number ; maximal external edge-number; total internal edge-weight; total external edge-weight; total internal edge-number ; total external edge-number)
        lst.append((Q, T, node, temp[0], temp[1], temp[2], temp[3]))
    '''sort by  Q decrease order and T increased order'''
    lst_add = sorted(lst, key=lambda i: (i[0], -int(i[1])), reverse=True)

    first_node = lst_add[0]
    number_min_temp = first_node[0]
    number_mout_temp = first_node[1]
    new_p = pv.p_upper_exact(n, N, M, first_node[0], first_node[1])

    if new_p < pvalue and abs(oldp - new_p) > gl.err_diff:
        add_node = first_node[2]
        pvalue = new_p
        wi = first_node[3]
        wo = first_node[4]
        ni = first_node[5]
        no = first_node[6]
        gl.item = gl.item + 1
    if first_node[1] == number_mout_minmum:
        return add_node, pvalue, wi, wo, ni, no
    length = len(lst_add)
    for i_index in range(1, length):
        i = lst_add[i_index]

        # if they have same T
        if i[1] == number_mout_minmum:
            new_p = pv.p_upper_exact(n, N, M, i[0], i[1])
            if new_p < pvalue and abs(oldp - new_p) > gl.err_diff:
                add_node = i[2]
                pvalue = new_p
                wi = i[3]
                wo = i[4]
                ni = i[5]
                no = i[6]
                gl.item = gl.item + 1
            break
        # it has more smaller T and unequal Q
        elif number_min_temp != i[0] and i[1] < number_mout_temp:
            new_p = pv.p_upper_exact(n, N, M, i[0], i[1])

            number_min_temp = i[0]
            number_mout_temp = i[1]

            if new_p < pvalue and abs(oldp - new_p) > gl.err_diff:
                pvalue = new_p
                add_node = i[2]
                wi = i[3]
                wo = i[4]
                ni = i[5]
                no = i[6]
                gl.item = gl.item + 1
    return add_node, pvalue, wi, wo, ni, no

def remove_function(pvalue_orig,ns,N,M,W_inc,W_de,S,edge_dict,W_ns_in,W_ns_out,number_ns_in,number_ns_out,last_pvalue):
    """
    for each node of ns, we try to find one node that removed from ns and make the p-value of the new ns is much lower than original p-value or the difference compared with addP not bigger than a threshold

    :param neib_dict: the neighbors dict (whose key is the node index and value is the number of neightbors in ns)
    :param pvalue: the p-value of the original ns
    :param ns: the candidate community node set
    :param N: total node number
    :param M: total edge number
    :param W_inc:  the sum list of the increased weight order . the increase order of the normalized edge weight is {w1, w2, ..., wM} where, w1 <= w2 <= ... <= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...+wM}
    :param W_de: similar to W_inc, the sum list of the decreased weight order . the decrease order of the normalized edge weight is {w1, w2, ..., wM} where, w1 >= w2 >= ... >= wM  thus W_inc is {w1, w1+w2, ... , w1+w2+...+wM}
    :param S: the node strength value dict, key is node index and value is node strength
    :param edge_dict:  e.g. {node1: {node2: {weight: value1}, {node3: {weight: value2 }}}}, node1 and node2 is connected with edge-weight value1, node1 and node3 is connnected with edge-weight value2
    :param W_ns_in: the total internal edge-weight
    :param W_ns_out: the total external edge-weight
    :param number_ns_in: the total internal edge number
    :param number_ns_out: the total external edge number
    :param last_pvalue: after add_function(), the minimal p-value
    :return: the node can be removed from ns; the p-value of new ns with the node removed from ns; wi(wo): is the internal(external) edge-weight of the new ns;  ni(no) is the internal(external) edge number of new ns
    """
    wi = W_ns_in
    wo = W_ns_out
    ni = number_ns_in
    no = number_ns_out
    remove_node = -1
    n = len(ns)-1

    if n ==-1:
        return -1,pvalue_orig,W_ns_in,W_ns_out,number_ns_in,number_ns_out
    node_information = {}
    '''if they have '''
    nin_nout = {}
    nin_win = {}
    nin_wout = {}
    # oldp = pvalue_orig
    pvalue = last_pvalue
    for i in ns:
        '''according to the original edge number, edge-weight values  and the ones connected with node i to adjuste the according value after add i to ns'''
        Wns_node, WN_ns_node, number_ns_node, number_N_ns_node = cal_node_Min_Mout(ns, i, edge_dict, S)
        iin = W_ns_in - Wns_node
        iout = W_ns_out - WN_ns_node + Wns_node
        number_iin = number_ns_in - number_ns_node
        number_iout = number_ns_out - number_N_ns_node + number_ns_node
        if number_iout < 0:
            exit()
            # Wns_node, WN_ns_node, number_ns_node, number_N_ns_node = cal_node_Min_Mout(ns, i, edge_dict, S)
        '''store the win, wout, niin, nout information of the new ns (after remove node i from ns) to node_information dict'''
        node_information[i] = (iin, iout, number_iin, number_iout)

        '''for two new candidate community a1(ns - node1), a2(ns - node2), if they have the same internal edge-number, we just need to calculate the one who has smaller external edge-number'''
        if number_iin in nin_nout:
            t1 = nin_nout[number_iin]
            if number_iout < t1[0]:
                nin_nout[number_iin] = [number_iout, i]
        else:
            nin_nout[number_iin] = [number_iout, i]
        '''for two new candidate community a1(ns - node1), a2(ns - node2), if they have the same internal edge-number, we just need to calculate the one who has bigger internal edge-weight'''
        if number_iin in nin_win:
            if iin > nin_win[number_iin][0]:
                nin_win[number_iin] = [iin, i]
        else:
            nin_win[number_iin] = [iin, i]
        '''for two new candidate community a1(ns - node1), a2(ns - node2), if they have the same internal edge-number, we just need to calculate the one who has smaller external edge-weight'''
        if number_iin in nin_wout:
            if iout < nin_wout[number_iin][0]:
                nin_wout[number_iin] = [iout, i]
        else:
            nin_wout[number_iin] = [iout, i]


    nodes_tobe_calculated = set()
    for nin in nin_nout:
        nodes_tobe_calculated.add(nin_nout[nin][1])
        nodes_tobe_calculated.add(nin_win[nin][1])
        nodes_tobe_calculated.add(nin_wout[nin][1])


    lst = []
    number_mout_minmum = float('inf')

    '''calculating the candidate nodes: if they have the same Q, we calculate the node who has the smallest T;
    if they have different Q,we calculate the node whose T is smaller than those Ts that have been touched
    until we meet the smallest T'''

    for node in nodes_tobe_calculated:
        temp = node_information[node]
        Qlow, gupper = pv.cal_Q_g(W_inc, W_de, temp[0], temp[1])

        Q = math.ceil(temp[2] * gl.alpha_for_edge + (min(Qlow, n * (n - 1) / 2)) * (1 - gl.alpha_for_edge))
        T = math.floor(temp[3] * gl.alpha_for_edge + gupper * (1 - gl.alpha_for_edge))
        number_mout_minmum = min(number_mout_minmum,T)
        # lst(minimual internal edge-number ; maximal external edge-number; total internal edge-weight; total external edge-weight; total internal edge-number ; total external edge-number)
        lst.append((Q,T,node,temp[0],temp[1],temp[2],temp[3]))
    '''sort by  Q decrease order and T increased order'''
    lst_remove = sorted(lst, key=lambda i: (i[0], -int(i[1])), reverse=True)


    # first node
    first_node = lst_remove[0]
    number_min_temp = first_node[0]
    number_mout_temp = first_node[1]
    new_p = pv.p_upper_exact(n, N, M, first_node[0], first_node[1])

    if new_p < pvalue or (last_pvalue - new_p) > - gl.err_diff:
        remove_node = first_node[2]
        pvalue = new_p
        wi = first_node[3]
        wo = first_node[4]
        ni = first_node[5]
        no = first_node[6]
        gl.item = gl.item + 1
    if first_node[1] == number_mout_minmum:
        return remove_node,pvalue,wi,wo,ni,no
    length = len(lst_remove)

    for i_index in range(1,length):
        i = lst_remove[i_index]
        # if they have same T
        if i[1] == number_mout_minmum:
            new_p = pv.p_upper_exact(n,N,M,i[0], i[1])

            if new_p < pvalue or (last_pvalue - new_p) > - gl.err_diff:
                remove_node = i[2]
                pvalue = new_p
                wi = i[3]
                wo = i[4]
                ni = i[5]
                no = i[6]
                gl.item = gl.item + 1
            break
        # it has more smaller T and unequal Q
        elif number_min_temp != i[0] and i[1] < number_mout_temp:
            new_p = pv.p_upper_exact(n,N,M,i[0], i[1])

            number_min_temp = i[0]
            number_mout_temp = i[1]

            if new_p < pvalue or (last_pvalue - new_p) > - gl.err_diff:
                pvalue = new_p
                remove_node = i[2]
                wi = i[3]
                wo = i[4]
                ni = i[5]
                no = i[6]
                gl.item = gl.item + 1
    return remove_node, pvalue, wi, wo, ni, no

def cal_node_Min_Mout(ns,node,edge_dict,S):
    """
    the total internal edge-weight connected with node (Wns_node); the total external edge-weight connected with node (WN_ns_node);
    to get the total internal edge number connected with node(number_ns_node); to get the total external edge number connected with node (number_N_ns_node_dir)

    :param ns:
    :param node:
    :param edge_dict:
    :param S:
    :return:
    """
    Wns_node = 0
    number_ns_node = 0
    nei = edge_dict[node]
    for i in nei:
        if i in ns:
            Wns_node += float(nei[i]['weight'])
            number_ns_node += 1
    WN_ns_node_dir = S[node] - Wns_node
    number_N_ns_node_dir = len(nei) - number_ns_node
    return Wns_node,WN_ns_node_dir,number_ns_node,number_N_ns_node_dir

def add_search_next(add_node,edge_dict,neib_dict,ns):
    """
    add the add_node to ns and update the neighbor dict

    :param add_node:
    :param edge_dict:
    :param neib_dict:
    :param ns:
    :return:
    """
    neib_dict.pop(add_node)
    nei = edge_dict[add_node]
    for i in nei:
        if i not in ns:
            if i in neib_dict:
                neib_dict[i] = neib_dict[i] + 1
            else:
                neib_dict[i] = 1
    ns.add(add_node)

def remove_search_next(remove_node,edge_dict,neib_dict,ns):
    """
    remove the remove_node from ns and update the neighbor dict

    :param remove_node:
    :param edge_dict:
    :param neib_dict:
    :param ns:
    :return:
    """
    nei = edge_dict[remove_node]
    nei_in_ns = 0
    for i in nei:
        if i in ns:
            nei_in_ns += 1
        else:
            if neib_dict[i] > 1:
                neib_dict[i] = neib_dict[i] - 1
            else:
                neib_dict.pop(i)
    neib_dict[remove_node] = nei_in_ns
    ns.remove(remove_node)


