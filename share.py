def get_neighborhood_weights_edgesnumber(ns,edge_dict):
    """
    :param ns: this is the candidate community set
    :param edge_dict: e.g. {node1: {node2: {weight: value1}, {node3: {weight: value2 }}}}, node1 and node2 is connected with edge-weight value1, node1 and node3 is connnected with edge-weight value2
    :return: the neighbors dict (whose key is the node index and value is the number of neightbors in ns); Win (Wout) is the total edge-weight in (external) ns; number_in (number_out) is the total number of edges in (external) ns;
    """
    neighborhood_dict = {}
    Win = 0
    number_in = 0
    Wout = 0
    number_out = 0
    for i in ns:
        neii = edge_dict[i]
        for j in neii:
            if j in ns:
                Win += float(neii[j]['weight'])
                number_in += 1
            else:
                Wout += float(neii[j]['weight'])
                number_out += 1
                if (j not in neighborhood_dict):
                    neighborhood_dict[j] = 1
                else:
                    neighborhood_dict[j] = neighborhood_dict[j] + 1
    Win /= 2
    number_in /= 2
    return neighborhood_dict,Win, number_in, Wout, number_out