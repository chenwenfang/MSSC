from scipy.special import comb
import datetime
import networkx as nx
#from sklearn import metrics
import re
def cal_purity(real,result):

    N = 0
    for i in range(0,len(real)):
        real[i] = set(real[i])

    for i in range(0,len(result)):
        N += len(result[i])
        result[i] = set(result[i])
        
    #print(N)
    M = 0

    for i in result:
        c = 0
        for j in real:

            temp = len(i & j)
            if temp > c:
                c = temp
        M += c

    return (M*1.0)/(N*1.0)

'''write by weixiaoqi'''
def cal_f1_randIndex(real,result):
    statis = []
    total = {}
    number_of_result = {}

    for i in range(0,len(real)):
        total[i] = len(real[i])
    result_len = []
    result_num = []
    label = {}
    TP = TN = FP = FN = 0
    TP_FP = 0
    TN_FN = 0
    for i in range(0,len(result)):
        result_len.append(len(result[i]))
# assign each node with its community index
    for i in range(0,len(real)):
        for j in real[i]:
            if j not in label:
                label[j] = set()
                label[j].add(i)
            else:
                label[j].add(i)
# get the number of nodes for each community
    for i in range(0,len(result)):
        for j in result[i]:
            if j not in label:
                label[j] = set("#")
            for lab in label[j]:
                if lab not in number_of_result:
                    number_of_result[lab] = 1
                else:
                    number_of_result[lab] += 1

    #print(number_of_result)
    #print(total)
    #print(result_len)
    #print(label)
# get the list of result community as [l(0)+l(1)+l(2)+..+l(n-1), l(1)+l(2)+..+l(n-1), ...., l(n-3)+l(n-2)+l(n-1), l(n-2) + l(n-1), 0]
    ax = 0
    for i in range(len(result_len)-1,-1,-1):
        ax += result_len[i]
        if i == len(result_len)-1:
            result_num.append(0)
        else:
            result_num.append(ax)
    result_num.reverse()

    #print(result_num)
    #print(result_len)
#get the number of nodes that in the groundTruth index in the result community
    for i in result:
        temp = dict()
        for j in i:
            if j not in label:
                continue
            for lab in label[j]:
                if lab in temp:
                    temp[lab] += 1
                else:
                    temp[lab] = 1
        statis.append(temp)
    #print(statis)

#TP+FP is the sum of each combinal number of each result community
    for i in result:
        TP_FP += comb(len(i),2)

    #for i in statis:
    #    for k in i:
    #        TP += comb(i[k],2)
    for nodes in result:
        for i in range(len(nodes)):
            for j in range(i+1,len(nodes)):
                if len(label[nodes[i]] & label[nodes[j]]) != 0:
                    #print(nodes[i],nodes[j])
                    TP+=1



    FP = TP_FP - TP

    #print(TP_FP,TP,FP)
    for i in range(0,len(result_num)):
        if result_num[i] == 0:
            break
        TN_FN += (result_num[i]-result_len[i])*result_len[i]

    #print(TN_FN)

    #print(total)
    #print(statis)
    for i in range(len(result)):
        for node_i in result[i]:
            for j in range(i+1,len(result)):
                for node_j in result[j]:
                    if( len(label[node_i] & label[node_j]) != 0):
                        FN += 1
    '''
    temp_total = copy.deepcopy(number_of_result)
    for i in statis:
        #print(i)
        for k in i:
            if temp_total[k] == 0:
                continue
            #print(FN,i[k],(temp_total[k]-i[k]))
            FN += i[k] *(temp_total[k]-i[k])
            temp_total[k] -= i[k]
    '''
    TN = TN_FN - FN
    precision = (TP*1.0)/((TP+FP)*1.0)
    recall = (TP*1.0)/((TP+FN)*1.0)

    f1 = 2.0*precision*recall/(recall+precision)
    rand_index = ((TP+TN)*1.0)/((TP+TN+FP+FN)*1.0)
    return f1,rand_index,precision,recall
def run(real_path, result_path):
    # real_path = "karate//karate-community.txt"
    # result_path = "result//karate_result.txt"

    f = open(real_path)
    lines = f.readlines()
    real = []
    for line in lines:
        if line == " " or line == "\n":
            continue
        data = re.split("\t| ", line)
        temp = []
        for j in data:
            if j == "" or j == '' or j == "\n":
                continue
            temp.append(j.strip())
        real.append(temp)
    f.close()

    f = open(result_path)
    lines = f.readlines()
    result = []
    for line in lines:
        if line == " " or line == "\n":
            continue
        data = re.split("\t| ", line)
        temp = []
        for j in data:
            if j == "" or j == '' or j == "\n":
                continue
            temp.append(j.strip())
        result.append(temp)

    f.close()
    community = {}

    for i in range(0, len(real)):
        for j in real[i]:
            if j not in community:
                community[j] = i
            else:
                continue

    f1, rand_index, precision, recall = cal_f1_randIndex(real, result)
    purity = cal_purity(real, result)

    return f1, rand_index, precision, recall, purity


def get_groundTruth_info(real_path):
    f = open(real_path)
    lines = f.readlines()
    real = []
    for line in lines:
        if line == " " or line == "\n":
            continue
        data = re.split("\t| ", line)
        temp = []
        for j in data:
            if j == "" or j == '' or j == "\n":
                continue
            temp.append(j.strip())
        real.append(temp)
    f.close()

    total = {}
    for i in range(0,len(real)):
        total[i] = len(real[i])
    # assign each node with its community index
    label = {}
    for i in range(0, len(real)):
        for j in real[i]:
            if j not in label:
                label[j] = set()
                label[j].add(i)
            else:
                label[j].add(i)
    return real, label
def get_result_info(result_path):
    f = open(result_path)
    lines = f.readlines()
    result = []
    for line in lines:
        if line == " " or line == "\n":
            continue
        data = re.split("\t| ", line)
        temp = []
        for j in data:
            if j == "" or j == '' or j == "\n":
                continue
            temp.append(j.strip())
        result.append(temp)
    f.close()
    return result
def cal_f1_randIndex_with_groundTruth(label,result):
    statis = []
    number_of_result = {}
    result_len = []
    result_num = []
    TP = TN = FP = FN = 0
    TP_FP = 0
    TN_FN = 0
    for i in range(0,len(result)):
        result_len.append(len(result[i]))

# get the number of nodes for each community
    for i in range(0,len(result)):
        for j in result[i]:
            if j not in label:
                label[j] = set("#")
            for lab in label[j]:
                if lab not in number_of_result:
                    number_of_result[lab] = 1
                else:
                    number_of_result[lab] += 1
# get the list of result community as [l(0)+l(1)+l(2)+..+l(n-1), l(1)+l(2)+..+l(n-1), ...., l(n-3)+l(n-2)+l(n-1), l(n-2) + l(n-1), 0]
    ax = 0
    for i in range(len(result_len)-1,-1,-1):
        ax += result_len[i]
        if i == len(result_len)-1:
            result_num.append(0)
        else:
            result_num.append(ax)
    result_num.reverse()
#get the number of nodes that in the groundTruth index in the result community
    for i in result:
        temp = dict()
        for j in i:
            if j not in label:
                continue
            for lab in label[j]:
                if lab in temp:
                    temp[lab] += 1
                else:
                    temp[lab] = 1
        statis.append(temp)
#TP+FP is the sum of each combinal number of each result community
    for i in result:
        TP_FP += comb(len(i),2)

    for nodes in result:
        for i in range(len(nodes)):
            for j in range(i+1,len(nodes)):
                if len(label[nodes[i]].intersection(label[nodes[j]])) != 0:
                    #print(nodes[i],nodes[j])
                    TP+=1
    FP = TP_FP - TP

    for i in range(0,len(result_num)):
        if result_num[i] == 0:
            break
        TN_FN += (result_num[i]-result_len[i])*result_len[i]

    for i in range(len(result)):
        for node_i in result[i]:
            for j in range(i+1,len(result)):
                for node_j in result[j]:
                    if( len(label[node_i] & label[node_j]) != 0):
                        FN += 1
    '''
    temp_total = copy.deepcopy(number_of_result)
    for i in statis:
        #print(i)
        for k in i:
            if temp_total[k] == 0:
                continue
            #print(FN,i[k],(temp_total[k]-i[k]))
            FN += i[k] *(temp_total[k]-i[k])
            temp_total[k] -= i[k]
    '''
    TN = TN_FN - FN
    precision = (TP * 1.0) / ((TP + FP) * 1.0)
    recall = (TP * 1.0) / ((TP + FN) * 1.0)

    f1 = 2.0 * precision * recall / (recall + precision)
    rand_index = ((TP + TN) * 1.0) / ((TP + TN + FP + FN) * 1.0)
    return f1,rand_index,precision,recall


# def cal_ARI(graphPath,realPath,resultPath):
#     G = nx.read_edgelist(graphPath, delimiter=' ', nodetype=int)
#     f = open(realPath)
#     lines = f.readlines()
#     real = []
#     for line in lines:
#         line = line.replace("\n","")
#         l = line.split(" ")
#         if "" in l:
#             l.remove("")
#         temp = []
#         for j in l:
#             temp.append(int(j))
#         real.append(temp)
#     f.close()
#
#     f = open(resultPath)
#     lines = f.readlines()
#     result = []
#     for line in lines:
#         if line == "":
#             continue
#         line = line.replace("\n","")
#         l = line.split(" ")
#         l.remove("")
#         if len(l) == 0:
#             continue
#         temp = []
#         for j in l:
#             temp.append(int(j))
#         result.append(temp)
#     f.close()
#     trues = {}
#     predict = {}
#     for i,j in enumerate(real):
#         for k in j:
#             trues[k] = i
#
#     for i,j in enumerate(result):
#         for k in j:
#             predict[k] = i
#     #print(predict)
#     #print(trues)
#     real_list = []
#     result_list = []
#     for i in G.nodes():
#         if i not in trues:
#             real_list.append(-1)
#         else:
#             real_list.append(trues[i])
#
#         if i not in predict:
#             result_list.append(-1)
#         else:
#             result_list.append(predict[i])
#     score = metrics.adjusted_rand_score(real_list,result_list)
#     print(score)
#     return score
def lfr():
    # LFR benchmark
    path1 = "./data/simulate_weighted_networks_realcommunity/"
    path2 = ["small_1000/", "big_1000/", "overlap/small_1000/", "overlap/big_1000/"]
    path3 = ["muw0.1/", "muw0.1/", "muw0.2/", "muw0.3/", "muw0.4/", "muw0.5/"]
    path4 = ["surprise.txt", "louvain.txt", "chinesewhispers.txt", "ipca.txt", "pycombo.txt", "rber_pots.txt",
             "mssc.txt", "infomap.txt", "slpaw.txt", "OSLOM.txt", "CCME_result.txt"]
    for i4 in path4:
        print(i4)
        for i2 in path2:
            for i3 in path3:
                truthcommunity = path1 + i2 + i3 + "community_handled.txt "
                resultCommunity = path1 + i2 + i3 + i4
                # print(resultCommunity)
                f1, rand_index, precision, recall, purity = run(truthcommunity, resultCommunity)
                # print(i2,i3, rand_index, f1, precision, recall, purity)
                print(i2, i3, purity)
def networkRepo():
    # networkRepo
    data1 = ["BZR-MD", "COX2-MD", "DHFR-MD", "ER-MD"]
    path4 = ["surprise.txt", "louvain.txt", "chinesewhispers.txt", "ipca.txt", "pycombo.txt", "rber_pots.txt",
             "mssc.txt", "infomap.txt", "slpaw.txt", "OSLOM.txt", "CCME.txt"]
    for i4 in path4:
        for i1 in data1:
            resultCommunity = "./data/real_data_sets/networkrepository/" + i1 + "_" + i4
            truthcommunity = "./data/real_data_sets/networkrepository/" + i1 + "_groundTruth.txt"
            f1, rand_index, precision, recall, purity = run(truthcommunity, resultCommunity)
            print(i4, i1, "f1, rand_index, precision, recall, purity: ", f1, rand_index, precision, recall, purity)
def weightedPPI():
    data1 = ["collins2007", "gavin2006_socioaffinities_rescaled", "krogan2006_core", "krogan2006_extended"]
    path4 = ["surprise.txt", "louvain.txt", "chinesewhispers.txt", "ipca.txt", "pycombo.txt", "rber_pots.txt",
            "mssc.txt","infomap.txt", "slpaw.txt", "oslom.txt", "ccme.txt"]
    path2 = ["CYC2008_complex.txt", "mips_3_100.txt", "sgd.txt"]
    for i2 in path2:
        for i4 in path4:
            for i1 in data1:
                resultCommunity = "./data/bio_datasets/index_term/" + i1 + "/" + i4
                truthcommunity = "./data/gold_standard/index_term/" + i1 + "/" + i2
                f1, rand_index, precision, recall, purity = run(truthcommunity, resultCommunity)
                print(i2, i4, i1, "rand_index,  purity: ",  rand_index,  purity)
def weightedUCI():
    # networkRepo
    data1 = ["iris", "ecoli", "wdbc", "movement_libras","newthyroid","wine","glass","appendicitis","monk-2"]
    path4 = ["surprise.txt", "louvain.txt", "chinesewhispers.txt", "ipca.txt", "pycombo.txt", "rber_pots.txt",
             "mssc.txt", "infomap.txt", "slpaw.txt", "oslom.txt", "ccme.txt"]
    for i4 in path4:
        for i1 in data1:
            if i4 == "ipca.txt" and (i1 == "wdbc" or i1 == "wine"):
                continue
            resultCommunity = "./data/weighted_real/UCI/" + i1 + "_" + i4
            truthcommunity = "./data/weighted_real/UCI/" + i1 + "_truth_community.txt"
            f1, rand_index, precision, recall, purity = run(truthcommunity, resultCommunity)
            print(i4, i1, "f1, rand_index, precision, recall, purity: ", f1, rand_index, precision, recall, purity)
def unweighted():
    # networkRepo
    data1 = ["karate", "football", "personal", "polblogs", "polbooks", "railway"]
    path4 = ["surprise.txt", "louvain.txt", "chinesewhispers.txt", "ipca.txt", "pycombo.txt", "rber_pots.txt",
             "mssc.txt", "infomap.txt", "slpaw.txt", "oslom.txt", "ccme.txt"]
    for i4 in path4:
        for i1 in data1:
            if i4 == "ccme.txt" and (i1 == "karate" or i1 == "football" ):
                continue
            resultCommunity = "./data/real_data_sets/" + i1 + "/" + i4
            truthcommunity = "./data/real_data_sets/" + i1 + "/" + i1 + "-community.txt"
            f1, rand_index, precision, recall, purity = run(truthcommunity, resultCommunity)
            print(i4, i1, "f1, rand_index, precision, recall, purity: ", f1, rand_index, precision, recall, purity)
def unweighted_bigdata():
    # networkRepo
    data1 = [ "dblp"]#["amazon", "dblp"]
    path4 = ["louvain.txt","surprise.txt",  "chinesewhispers.txt", "rber_pots.txt", "infomap.txt", "slpaw.txt", "oslom.txt","mssc.txt"]#, "ccme.txt", "ipca.txt", "pycombo.txt"]
    for i1 in data1:
        truthcommunity = "./data/real_data_sets/" + i1 + "/" + i1 + "-community.txt"
        real, label = get_groundTruth_info(truthcommunity)
        for i4 in path4:
            print(i4,"start: ",datetime.datetime.now())
            resultCommunity = "./data/real_data_sets/" + i1 + "/" + i4
            result = get_result_info(resultCommunity)
            f1, rand_index, precision, recall = cal_f1_randIndex_with_groundTruth(label, result)
            purity = cal_purity(real, result)
            print(i4,  "f1, rand_index, precision, recall, purity: ", f1, rand_index, precision, recall, purity)
            print(i4,"end Time:", datetime.datetime.now())
def main():
    unweighted_bigdata()
    # lfr()
    # networkRepo()
    # unweighted()
    # weightedUCI()
    # weightedPPI()


if __name__ == '__main__':
    main()
    # run("real.txt","result.txt")


'''write by ChenWenFang'''
# this program needs more memerory
# def cal_f1_randIndex(real, result):
#     tp_fp = 0
#     # get the node pairs of the whole groundTruth communities
#     node_pairs = set()
#     for ireal in real:
#         len_ith_real = len(ireal)
#         for node1 in range(0,len_ith_real):
#             for node2 in range(node1+1, len_ith_real):
#                 node_pairs.add((ireal[node1],ireal[node2]))
#
#     tp = 0
#     for index_result in range(0,len(result)):
#         iresult = result.get(index_result)
#         for noderes1 in range(0,len(iresult)):
#             tp_fp = tp_fp + len(iresult.get(noderes1))*(len(iresult.get(noderes1))-1)/2
#             for noderes2 in range(noderes1+1, len(iresult)):
#                 if (iresult[noderes1],iresult[noderes2]) in node_pairs or (iresult[noderes2],iresult[noderes1]) in node_pairs:
#                     tp = tp + 1
#     tn_fn = 0
#     tn = 0
#     result_node_pairs = []
#     for index_result in range(0,len(result)):
#         iresult = result.get(index_result)
#         for index_result2 in range(index_result,len(result)):
#             tn_fn = tn_fn + len(iresult)*len(result.get(index_result2))
#             iresult2 = result[index_result2]
#             for inode1 in iresult:
#                 for inode2 in iresult2:
#                     if inode1 != inode2 and (inode1,inode2) not in node_pairs and (inode2,inode1) not in node_pairs:
#                         tn = tn + 1
#     fn = tn_fn - tn
#     precision = tp / tp_fp
#     recall = tp / (tp + fn)
#
#     f1 = 2 * precision * recall / (recall + precision)
#     rand_index = (tp + tn) / (tn_fn + tp_fp)
#     return precision,recall,f1,rand_index
