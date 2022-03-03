import networkx as nx
import global_var_model as gl

#-------------------normalize the weights
def init_simulate_graph_txt(file,normalize_type,weighted_type):
    """
    :param file: the input data file with three columns, the first one is the node , the second one is the node, the third one is the edge-weight
    :param normalize_type: the way of normalizing original-edge-weight
    :return:
    """
    fp = open(file,"r");lines = fp.readlines(); fp.close()
    index = 1
    #---------------- nodes keep the original names of all vertexs
    nodes = set()
    global net_old
    net_old = nx.Graph()
    global net
    net = nx.Graph()
    node_node_weight = [];new_node_node_weight = []
    coefficient_degree = []
    degrees = []
    coefficient = []
    S = {}; S_old = {}
    Wlist = []; Wlist_old = [];
    link_table = [];
    link_table.append([])
    # ----------------- read data files and get the name-index dict, index-name dict and the original network
    # if the data files have the weight column, we use the weight value, else, we set the original weight is 1.0
    if(weighted_type == False):
        gl.alpha_for_edge = 1
        normalize_type = "NO"
        for line in lines:
            data = line.split()
            nodestr1 = data[0].strip()
            nodestr2 = data[1].strip()
            weight = 1.0
            if nodestr1 != nodestr2:
                if nodestr1 not in gl.name_index:
                    gl.name_index[nodestr1] = index
                    gl.index_name[index] = nodestr1
                    index += 1
                    link_table.append([])
                node1 = gl.name_index[nodestr1]
                if nodestr2 not in gl.name_index:
                    gl.name_index[nodestr2] = index
                    gl.index_name[index] = nodestr2
                    index += 1
                    link_table.append([])
                node2 = gl.name_index[nodestr2]
                nodes.add(nodestr1)
                nodes.add(nodestr2)
                node_node_weight.append((node1,node2,weight))
    else:
        for line in lines:
            data = line.split()
            nodestr1 = data[0].strip()
            nodestr2 = data[1].strip()
            weight = float(data[2].strip())
            if nodestr1 != nodestr2:
                if nodestr1 not in gl.name_index:
                    gl.name_index[nodestr1] = index
                    gl.index_name[index] = nodestr1
                    index += 1
                    link_table.append([])
                node1 = gl.name_index[nodestr1]
                if nodestr2 not in gl.name_index:
                    gl.name_index[nodestr2] = index
                    gl.index_name[index] = nodestr2
                    index += 1
                    link_table.append([])
                node2 = gl.name_index[nodestr2]
                nodes.add(nodestr1)
                nodes.add(nodestr2)
                node_node_weight.append((node1, node2, weight))

    net_old.add_weighted_edges_from(node_node_weight)
    pre_edge_list_old = nx.to_dict_of_dicts(net_old)

    if normalize_type == "simple":
        #---------------- getting the original strength of each node
        for id,val in pre_edge_list_old.iteritems():
            strength = 0.0
            for j in val:
                strength += float(val[j]['weight'])
            S_old[id] = strength
        #----------------- getting the normalized weight, degree[stored in terms of list((node,degree))] and strength[stored in terms of dict with node is the key] of each edge
        #------------------ normalized_weight of node i and j = (original weight of node i and j) / (original strength of node i + original strength of node j)
        for id,val in pre_edge_list_old.iteritems():
            new_strength = 0.0
            leng = len(val)
            degrees.append((id, leng))
            for j in val:
                link_table[id].append(j)
                weight_old = val[j]['weight']
                weight = float(weight_old) / (S_old[id] + S_old[j])
                new_node_node_weight.append((id, j, weight))
                new_strength += weight
                if j > id:
                    Wlist.append(weight)
            S[id] = new_strength
    elif normalize_type == "NO":
        for id, val in pre_edge_list_old.iteritems():
            strength = 0.0
            for j in val:
                link_table[id].append(j)
                enheng = float(val[j]['weight'])
                strength += enheng
                if j > id:
                    Wlist.append(enheng)
            S[id] = strength
            leng = len(val)
            degrees.append((id, leng))
        new_node_node_weight = node_node_weight
    net.add_weighted_edges_from(new_node_node_weight)
    pre_edge_list = nx.to_dict_of_dicts(net)
    #-------------------- calculate the combination of coefficient and degree [coefficient fenzi /degree of node i] as the node selected standard
    for id in pre_edge_list_old:
        coeff_fenzi = 0
        name = gl.index_name[id]
        temp = link_table[id]
        leng = len(temp)
        if leng <= 3:
            coeffi_degree = 0
            coeffi = 0
        else:
            for j in range(leng):
                for k in range(j+1,leng):
                    if temp[k] in pre_edge_list_old[temp[j]]:
                        coeff_fenzi += 1
            coeffi_degree = coeff_fenzi/(leng*1.0)
            coeffi = coeff_fenzi/(leng*(leng-1.0)/2.0)
        coefficient.append((id,coeffi))
        coefficient_degree.append((id,coeffi_degree))
    return net,pre_edge_list,coefficient_degree,link_table,Wlist, S, nodes,degrees,coefficient


