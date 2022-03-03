import init_graph
import search
import share
import math
import copy
import global_var_model as gl
import pv
import sys

def get_groundTruth_resultCommunity_information():
    log_factorial()
    # path1 = "./real_datasets/networkrepository/"
    # path2 = ["FIRSTMM-DB"]#,"BZR-MD", "DHFR-MD", "ER-MD" , "COX2-MD","FIRSTMM-DB"]
    # for i2 in path2:
    #     path_input = path1 + i2 + "_countdown.txt"
    #     weighted_type = True
    #     path_groundTruth = path1 + i2 + "_groundTruth.txt"

    # path1 = ["real_datasets\\moody\\moody"]
    # for i2 in path1:
    #     path_input = i2 + "_commun.txt"
    #     path_output_lh = i2 + "_commun_MSSC.txt"
    #     weighted_type = True
    #     path_groundTruth = i2 + "_groundTruth.txt"

    #  # this is the input path for simulated overlap datasets with weight controlled networks whose muw ranges from 0 to 0.4 and mut fixed at 0.5
    # path1 = "simulate_weighted_networks_realcommunity\\overlap\\big_5000\\"
    # # path2 = ["big_5000\\"]  # ["small_1000","big_1000\\","small_5000","big_5000"]#,"big_5000"]
    # path4 = ["muw0mut0.5\\"]#, "\\muw0.1mut0.5\\","\\muw0.2mut0.5\\","\\muw0.3mut0.5\\","\\muw0.4mut0.5\\"]
    # path5 = "network.dat"
    # path7 = "MSSC.txt"
    # for i2 in path4:
    #     path_input = path1 + i2 + path5
    #     path_output_lh = path1 + i2 + path7
    #     weighted_type = True
    #     path_groundTruth = path1 + i2 + "community_handled.txt"

    # # #  # this is the input path for simulated overlap datasets
    path1 = "simulate_weighted_networks_realcommunity\\overlap\\big_1000\\"
    # path1_result = "result_simulate_weighted_networks_realcommunity\\overlap\\small_1000\\"
    # path2 = ["small_1000\\"]#["small_1000","big_1000\\","small_5000","big_5000"]#,"big_5000"]
    path4 = ["muw0\\","\\muw0.1\\","\\muw0.2\\","\\muw0.3\\","\\muw0.4\\","\\muw0.5\\"]
    path5 = "network.dat"
    path7 = "lh_remove_result.txt"
    for i2 in path4:
        path_input = path1 + i2 + path5
        # path_output_lh = path1 + i2 + path7
        weighted_type = True
        path_groundTruth = path1 + i2 + "community_handled.txt"

    #  # this is the input path for simulated unoverlap datasets
    # path1 = "simulate_weighted_networks_realcommunity\\small_1000\\"
    # path1_result = "result_simulate_weighted_networks_realcommunity\\small_1000\\"
    # # path2 = ["small_1000"]#, "big_1000\\", "small_5000", "big_5000"]  # ,"big_5000"]
    # path4 = ["\\muw0.4\\"]#, "\\muw0.1\\", "\\muw0.2\\", "\\muw0.3\\", "\\muw0.4\\", "\\muw0.5\\"]
    # path5 = "network.dat"
    # path7 = "lh_remove_result.txt"
    # for i2 in path4:
    #     path_input = path1 + i2 + path5
    #     path_output_lh = path1 + i2 + path7
    #     weighted_type = True
    #     path_groundTruth = path1 + i2 + "community_handled.txt"

        gl.index_name.clear()
        gl.name_index.clear()
        net, edge_dict, coefficient_degree, link_table, Wlist, S, nodes, degrees, coefficient = init_graph.init_simulate_graph_txt(path_input, "NO", weighted_type)
        Wlist_de = sorted(Wlist, reverse=True)
        Wlist_inc = sorted(Wlist, reverse=False)
        Wlist_decrease = []
        Wlist_increase = []
        temp_de = 0
        temp_inc = 0
        M = net.number_of_edges()
        N = net.number_of_nodes()
        # print M/(N * (N - 1.0) / 2.0), (2.0) * M / (N * 1.0)

        for i in range(M):
            temp_de += Wlist_de[i]
            Wlist_decrease.append(temp_de)
            temp_inc += Wlist_inc[i]
            Wlist_increase.append(temp_inc)

        '''get the information of the result communities'''
        path_output_lh = path1 + i2 + "MSSC.txt"
        file_result = open(path_output_lh,'r')
        lines = file_result.readlines()
        print len(gl.name_index)
        for line in lines:
            data = line.split()
            curCom = []
            # print data
            for idata in data:
                curCom.append(gl.name_index[(idata.strip())])
            ns = set(curCom)
            nb_dict, W_ns_in, number_ns_in, W_ns_out, number_ns_out = share.get_neighborhood_weights_edgesnumber(ns,edge_dict)
            '''get the pvalue of the community'''
            n = len(ns)
            Qlow, gupper = pv.cal_Q_g(Wlist_increase,Wlist_decrease, W_ns_in, W_ns_out)
            '''according to the parameter of weight influence and edge-number influence adjust the exact edge-number'''
            T = math.floor(gl.alpha_for_edge * number_ns_out + (1 - gl.alpha_for_edge) * (gupper))
            Q = math.ceil(gl.alpha_for_edge * number_ns_in + (1 - gl.alpha_for_edge) * (min(Qlow, n * (n - 1) / 2)))

            pvalue = pv.p_upper_exact(n, N, M, Q, T)
            # print "result commuity",sorted(ns)
            a1 = n * (n - 1.0) / 2.0
            a2 = n * (N - n) * 1.0
            # print "groundTruth commuity", sorted(nsg)
            print "W_ns_in, number_ns_in, W_ns_out, number_ns_out,pvalue, av_in_we, av_ex_we, pin, pout", W_ns_in, number_ns_in, W_ns_out, number_ns_out, pvalue, W_ns_in / number_ns_in, W_ns_out / number_ns_out, number_ns_in / a1, number_ns_out / a2

        '''get the information of the groundTruth communities'''
        file_groundTruth = open(path_groundTruth, 'r')
        linesg = file_groundTruth.readlines()
        for lineg in linesg:
            datag = lineg.split()
            curComg = []
            for idatag in datag:
                curComg.append(gl.name_index[idatag.strip()])
            nsg = set(curComg)
            nb_dictg, W_ns_ing, number_ns_ing, W_ns_outg, number_ns_outg = share.get_neighborhood_weights_edgesnumber(nsg,
                                                                                                                 edge_dict)
            '''get the pvalue of the community'''
            ng = len(nsg)
            Qlowg, gupperg = pv.cal_Q_g(Wlist_increase, Wlist_decrease, W_ns_ing, W_ns_outg)
            '''according to the parameter of weight influence and edge-number influence adjust the exact edge-number'''
            Tg = math.floor(gl.alpha_for_edge * number_ns_outg + (1 - gl.alpha_for_edge) * (gupperg))
            Qg = math.ceil(gl.alpha_for_edge * number_ns_ing + (1 - gl.alpha_for_edge) * (min(Qlowg, ng * (ng - 1) / 2)))

            pvalueg = pv.p_upper_exact(ng, N, M, Qg, Tg)
            a11 = ng*(ng-1.0)/2.0
            a22 = ng*(N-ng)*1.0
            # print "groundTruth commuity", sorted(nsg)
            # print "W_ns_in, number_ns_in, W_ns_out, number_ns_out,pvalue, av_in_we, av_ex_we, pin, pout", W_ns_ing, number_ns_ing, W_ns_outg, number_ns_outg, pvalueg, W_ns_ing/number_ns_ing, W_ns_outg/number_ns_outg, number_ns_ing/a11, number_ns_outg/a22
            print "W_ns_in, number_ns_in, W_ns_out, number_ns_out,pvalue, av_in_we, av_ex_we, pin, pout", W_ns_ing, number_ns_ing, W_ns_outg, number_ns_outg, pvalueg,number_ns_ing/a11,number_ns_outg/a22,ng
        print M/(N*(N-1)/2)

# ----------------prestore the log(N2!) into log_factorial_num
def log_factorial():
    """
    pre store log(upperLog!)
    :return:
    """
    gl.factorial_num.append(0.0)
    for i in range(1, gl.maximal_store):
        gl.factorial_num.append(gl.factorial_num[i - 1] + math.log(i))
def get_log_factorial(a,b):
    """
    calculate  log(n!)
    :param n:
    :return:
    """
    if a == 0 or b == 0:
        return 0
    elif a == b:
        return 0
    elif a > 10000 or b > 10000:
        return (0.5 + a) * math.log(a) - (0.5 + b) * math.log(b) - (0.5 + a - b) * math.log(a - b) - 0.5 * math.log(
            2.0 * math.pi)
    else:
        return gl.factorial_num[int(a)] - gl.factorial_num[int(b)] - gl.factorial_num[int(a - b)]

def run_MSSC(path_input,path_output,weighted_type,node_type):
    log_factorial()
    gl.name_index.clear()
    gl.index_name.clear()
    print ("input", path_input)
    without_remove = []
    '''
    net is the graph
    edge_dict is the node1{link_node1{normalized weight1},link_node2{normalized weight2}...}
    pre_edge_list_old is the node1{link_node1{original weight1},link_node2{original weight2}...}
    coefficient_degree is the seed node select standard in terms of dict, key is the node index
    link_table is the neighbor list of each node recored in terms of index
    Wlist is the list of all normalized weights if the data is normailzed
    S is the set of strength S[i] = strength of node i
    nodes is the set of the original names of all nodes;
    name_index={nodename,nodeindex}; index_name = {nodeindex,nodename}
    degree is the list((nodeindex,degree))
    get the formated network based on the input data
    '''
    net,edge_dict,coefficient_degree,link_table, Wlist, S, nodes,degrees,coefficient = init_graph.init_simulate_graph_txt(path_input,"simple",weighted_type)
    N = net.number_of_nodes()
    M = len(Wlist)
    average_degree = (2.0) * M / (N * 1.0)

    #------------choose the seed node selection standard
    if (average_degree >= 10) or (average_degree > 5 and average_degree < 10 and N > 1000):
        choose_bar = coefficient_degree
        gl.err_diff = 7
    else:
        choose_bar = degrees
        gl.err_diff = 0
   
    choose_bar.sort(key=lambda i: i[1], reverse=True)
    #------------choose_bar_dict is for finding the node is removed or not
    choose_bar_dict = {}
    for i in choose_bar:
        choose_bar_dict[i[0]] = i[1]
    #-----------getting the seed node and candidate community
    def choose_node(choose_bar):
        # get the node index who has the maximal coefficient or coefficient_degree value from the left nodes
        max_pos = choose_bar[0][0]
        # remove the node from the left nodes
        choose_bar.remove((max_pos, choose_bar[0][1]))
        # mark it in the dict form
        choose_bar_dict[max_pos] = -1
        # get the candidate community
        nslist = copy.deepcopy(link_table[max_pos])
        nslist.append(max_pos)
        ns = set(nslist)
        return ns, max_pos

    #-----------sort Wlist and store the sum of first i elements into Wlist_decrease and Wlist_increase seperatly
    Wlist_de = sorted(Wlist,reverse=True)
    Wlist_inc = sorted(Wlist, reverse=False)
    Wlist_decrease = []
    Wlist_increase = []
    temp_de = 0
    temp_inc = 0
    for i in range(M):
        temp_de += Wlist_de[i]
        Wlist_decrease.append(temp_de)
        temp_inc += Wlist_inc[i]
        Wlist_increase.append(temp_inc)
    while len(choose_bar) > 2:
        ns,max_pos = choose_node(choose_bar)
        if len(ns) > 1:
            #-------------------considering the redudency process of cal_MinMout, I combine search.call_MinMout with share.choose_neighborhood process
            nb_dict,W_ns_in, number_ns_in, W_ns_out, number_ns_out = share.get_neighborhood_weights_edgesnumber(ns,edge_dict)
            gl.item = 0
            sc,pvalue = search.search_community(nb_dict,ns,W_ns_in, number_ns_in, W_ns_out, number_ns_out,edge_dict,link_table,N,M,Wlist_increase,Wlist_decrease,S)
            '''add result to wihtout_remove list and remove nodes of sc from the left_nodes dict'''
            if len(sc)>2 and pvalue < gl.log_alpha:
                without_remove.append(sc)
                for i in sc:
                    if choose_bar_dict[i] >= 0.0:
                        choose_bar.remove((i,choose_bar_dict[i]))
                        choose_bar_dict[i] = -1
    #------------------clusterone overlap rate and remove redundency
    with_lh_remove_result = []
    with_lh_remove_result_temp = copy.deepcopy(without_remove)
    i = 0
    while i < len(with_lh_remove_result_temp):
        combin_flag = 0
        j = 0
        while j < i:
            ti = with_lh_remove_result_temp[i]
            tj = with_lh_remove_result_temp[j]
            temp_i = set(ti)
            temp_j = set(tj)
            l1 = len(temp_i & temp_j) * 1.0
            l2 = len(temp_i)
            l3 = len(temp_j)
            cover_liang = l1/min(l2,l3)
            if cover_liang >= gl.cover_rate:
                temp = temp_i | temp_j
                with_lh_remove_result_temp.remove(ti)
                with_lh_remove_result_temp.remove(tj)
                with_lh_remove_result_temp.append(list(temp))
                combin_flag = 1
                break
            j = j + 1
        if combin_flag == 1:
            i = i - 1
        else:
            i = i + 1
    for i in with_lh_remove_result_temp:
        if not with_lh_remove_result.__contains__(i):
            with_lh_remove_result.append(i)

    removed_result = []
    '''transform the index value to name'''
    leng = len(with_lh_remove_result)
    for i in range(leng):
        lala = with_lh_remove_result[i]
        sc_name = []
        for j in lala:
            sc_name.append(gl.index_name[j])
        removed_result.append(sc_name)


    f = open(path_output,'w')
    for i in removed_result:
        '''if the node name is int type, we can sort them in order'''
        if node_type == "int":
            ii = list(map(int,i))
            iii = sorted(ii)
            '''if the node name is a string, we just write them '''
        else:
            iii = list(i)
        leng4 = len(iii)
        for j in range(leng4):
            f.write(str(iii[j]))
            f.write(' ')
        f.write('\n')
    f.close()
    print "end"

def main(argv):
    path_input = argv[1]
    path_output = path_input + "_mssc.txt"
    weighted_type = argv[2]
    node_type = argv[3]
    print("netwrok",path_input)
    print("weighted_type",weighted_type)
    print("node_type",node_type)
    run_MSSC(path_input,path_output,weighted_type,node_type)

if __name__ == '__main__':
    main(sys.argv)
