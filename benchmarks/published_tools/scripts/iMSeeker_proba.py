#!/usr/bin/env python

import os
import argparse
import regex
import numpy as np
import pickle
from imblearn.ensemble import BalancedRandomForestClassifier
import xgboost



def input_chr_file(file):
    chr_dic={}
    with open(file,"r") as f:
        sig=0
        for line in f:
            seq=line.strip()
            if len(seq)==0:
                continue
            if seq[0]==">":
                if sig==1:
                    chr_dic[char_name]=Seq_total.upper().replace('U', 'T')
                char_name=seq[1:]
                Seq_total=""
                sig=1
            else:
                Seq_total=Seq_total+seq
        chr_dic[char_name]=Seq_total.upper().replace('U', 'T')
    return chr_dic


def merge_sort(lst,loc):
    if len(lst)<= 1:
        return lst
    mid = len(lst)//2
    left = merge_sort(lst[:mid],loc) 
    right = merge_sort(lst[mid:],loc)
    merged = [] 
    while left and right:
        if left[0][loc] <= right[0][loc]:
            merged.append(left.pop(0))
        else:
            merged.append(right.pop(0))
    merged.extend(right if right else left) 
    return merged


def combined_iMregion_and_NoniMregion(iMregion,Gene_range):
    sorted_iMregion=merge_sort(iMregion,1)
    iMregion_new=[]
    NoniMregion=[]
    for i in sorted_iMregion:
        if iMregion_new==[]:
            iMregion_new.append(i)
        else:
            not_need=[]
            test=i
            for j in iMregion_new:
                if j[0]<=test[0] and j[1]>=test[0]:
                    test=[min([j[0],test[0]]),max([j[1],test[1]])]
                    not_need.append(j)
                elif j[0]<=test[1] and j[1]>=test[1]:
                    test=[min([j[0],test[0]]),max([j[1],test[1]])]
                    not_need.append(j)
                elif j[0]>=test[0] and j[1]<=test[1]:
                    test=[min([j[0],test[0]]),max([j[1],test[1]])]
                    not_need.append(j)
                elif j[0]<=test[0] and j[1]>=test[1]:
                    test=[min([j[0],test[0]]),max([j[1],test[1]])]
                    not_need.append(j)
            iMregion_new.append(test)
            for j in not_need:
                iMregion_new.remove(j)
    iMregion_new=merge_sort(iMregion_new,1)
    if iMregion_new[0][0]-Gene_range[0]>=1:
        NoniMregion.append([Gene_range[0],iMregion_new[0][0]])
    if Gene_range[-1]-iMregion_new[-1][-1]>=1:
        NoniMregion.append([iMregion_new[-1][-1],Gene_range[-1]])
    for i in range(len(iMregion_new)-1):
        if iMregion_new[i+1][0]-iMregion_new[i][1]>=1:
            NoniMregion.append([iMregion_new[i][1],iMregion_new[i+1][0]])
    NoniMregion=merge_sort(NoniMregion,1)
    return iMregion_new, NoniMregion

#checked
def check_iM(piM):
    if piM.count("A")+piM.count("G")+piM.count("C")+piM.count("T")+piM.count("U")==len(piM):
        return True
    return False

#checked
def find_all_motifs(start_node,edge_dic,lis,result):
    if start_node in edge_dic.keys():
        abandon = []
        for j in edge_dic[start_node].keys():
            new=[]
            for n in lis:
                if n[-1]==start_node:
                    new.append(n+[j])
                    abandon.append(n)
            for n in new:
                if len(n)==4:
                    value=int(n[-1].split("-")[-1])-int(n[0].split("-")[0])
                    if n in lis:
                        abandon.append(n)
                    if n not in result:
                        result.append(n)
                else:
                    if n not in lis:
                        lis.append(n)
        for n in abandon:
            if n in lis:
                lis.remove(n)
        next_l=[]
        for n in lis:
            next_l.append(n[-1])
        for j in next_l:
            find_all_motifs(j, edge_dic, lis,result)




def identifier_imotif_directedGraphic_noLoopConstrain_alliM(seq, n, GCtype):
    C_track_tmp = "["+GCtype+"]{" + str(n) + "}"
    res = regex.finditer(C_track_tmp, seq, overlapped=True)
    c_list=[]
    for j in res:
        c_list.append([j.span()[0], j.span()[1]])
    c_dic={}
    for i in range(len(c_list)):
        c_dic[str(c_list[i][0])+"-"+str(c_list[i][1])]=i
    c_dic_re={}
    for i in range(len(c_list)):
        c_dic_re[str(i)]=str(c_list[i][0])+"-"+str(c_list[i][1])
    edge={}
    node=[]
    for i in range(len(c_list)):
        for j in range(len(c_list)):
            if c_list[i][0]-c_list[j][1]>=0:
                start=str(c_list[j][0])+"-"+str(c_list[j][1])
                end=str(c_list[i][0])+"-"+str(c_list[i][1])
                if start not in edge.keys():
                    edge[start]={}
                edge[start][end]=c_list[i][0]-c_list[j][1]
                if start not in node:
                    node.append(start)
    result_all=[]
    for i in node:
        result=[]
        lis=[[i]]
        find_all_motifs(i, edge, lis, result)
        result_all=result_all+result
    return result_all


def identifier_imotif_directedGraphic_alliM(seq, n, GCtype,loop_lower_bound,loop_higher_bound):
    C_track_tmp = "["+GCtype+"]{" + str(n) + "}"
    res = regex.finditer(C_track_tmp, seq, overlapped=True)
    c_list=[]
    for j in res:
        c_list.append([j.span()[0], j.span()[1]])
    c_dic={}
    for i in range(len(c_list)):
        c_dic[str(c_list[i][0])+"-"+str(c_list[i][1])]=i
    c_dic_re={}
    for i in range(len(c_list)):
        c_dic_re[str(i)]=str(c_list[i][0])+"-"+str(c_list[i][1])
    edge={}
    node=[]
    for i in range(len(c_list)):
        for j in range(len(c_list)):
            if c_list[i][0]-c_list[j][1]>=loop_lower_bound and c_list[i][0]-c_list[j][1]<=loop_higher_bound:
                start=str(c_list[j][0])+"-"+str(c_list[j][1])
                end=str(c_list[i][0])+"-"+str(c_list[i][1])
                if start not in edge.keys():
                    edge[start]={}
                edge[start][end]=c_list[i][0]-c_list[j][1]
                if start not in node:
                    node.append(start)
    result_all=[]
    for i in node:
        result=[]
        lis=[[i]]
        find_all_motifs(i, edge, lis, result)
        result_all=result_all+result
    return result_all


def choose_long_stem(select_tmp):
    if len(select_tmp)==1:
        return select_tmp
    max_stem=0
    all_dic={}
    for i in select_tmp:
        len_tmp=int(i[0].split("-")[1])-int(i[0].split("-")[0])
        if max_stem<len_tmp:
            max_stem=len_tmp
        if len_tmp not in all_dic.keys():
            all_dic[len_tmp]=[]
        all_dic[len_tmp].append(i)
    
    #for i in all_dic.keys():
        #print("choose_long_stem",i, all_dic[i])
    #print(max_stem,all_dic[max_stem])
    
    return all_dic[max_stem]


def choose_longest_conformation(select_tmp):
    select_tmp_long_stem = choose_long_stem(select_tmp)
    if len(select_tmp_long_stem)==1:
        return select_tmp_long_stem
    long=0
    all_dic={}
    for i in select_tmp_long_stem:
        len_tmp=int(i[-1].split("-")[-1])-int(i[0].split("-")[0])
        if long<len_tmp:
            long=len_tmp
        if len_tmp not in all_dic.keys():
            all_dic[len_tmp]=[]
        all_dic[len_tmp].append(i)
    return all_dic[long]

"""
def Find_longest_stem_for_choose_shortest_conformation(select_tmp):
    max_stem=0
    all_dic={}
    for i in select_tmp:
        len_tmp=int(i[0].split("-")[1])-int(i[0].split("-")[0])
        if max_stem<len_tmp:
            max_stem=len_tmp
        if len_tmp not in all_dic.keys():
            all_dic[len_tmp]=[]
        all_dic[len_tmp].append(i)
    return max_stem

def Detect_upper_boundary_for_choose_shortest_conformation(select_dic,short):
    cand=list(select_dic.keys())
    upp=max(cand)
    for i in range(0,upp+1):
        if short+i not in cand:
            return i-1
        
def choose_shortest_conformation(select_tmp):
    if len(select_tmp)==1:
        return select_tmp
    short=1000000000000000000000
    all_dic={}
    for i in select_tmp:
        len_tmp=int(i[-1].split("-")[-1])-int(i[0].split("-")[0])
        if short>len_tmp:
            short=len_tmp
        if len_tmp not in all_dic.keys():
            all_dic[len_tmp]=[]
        all_dic[len_tmp].append(i)
        
    upp=Detect_upper_boundary_for_choose_shortest_conformation(all_dic,short)
    if upp==0:
        return all_dic[short]
    else:
        stem_short=Find_longest_stem_for_choose_shortest_conformation(all_dic[short])
        #print(stem_short)
        for i in range(1,upp+1):
            stem_tmp=Find_longest_stem_for_choose_shortest_conformation(all_dic[short+i])
            if stem_short<stem_tmp:
                stem_short=stem_tmp
                final=i
            else:
                final=i-1
                break
        return all_dic[short+final]
"""
"""
def choose_shortest_conformation(select_tmp):
    if len(select_tmp)==1:
        return select_tmp
    short=1000000000000000000000
    all_dic={}
    for i in select_tmp:
        len_tmp=int(i[-1].split("-")[-1])-int(i[0].split("-")[0])-4*(int(i[0].split("-")[1])-int(i[0].split("-")[0]))
        if short>len_tmp:
            short=len_tmp
        if len_tmp not in all_dic.keys():
            all_dic[len_tmp]=[]
        all_dic[len_tmp].append(i)
    result=choose_long_stem(all_dic[short])
    return result
"""
def choose_shortest_conformation(select_tmp):
    select_tmp_long_stem = choose_long_stem(select_tmp)
    if len(select_tmp_long_stem)==1:
        return select_tmp_long_stem
    short=1000000000000000000000
    all_dic={}
    for i in select_tmp_long_stem:
        len_tmp=int(i[-1].split("-")[-1])-int(i[0].split("-")[0])
        if short>len_tmp:
            short=len_tmp
        if len_tmp not in all_dic.keys():
            all_dic[len_tmp]=[]
        all_dic[len_tmp].append(i)
    return all_dic[short]
        
def choose_min_std(select_tmp):
    if len(select_tmp)==1:
        #print("choose_min_std",select_tmp[0])
        return select_tmp[0]
    all_dic={}
    all_le=[]
    all_std=[]
    for i in select_tmp:
        tmp = []
        for j in range(len(i)-1):
            tmp.append(int(i[j+1].split("-")[0])-int(i[j].split("-")[1]))
        all_le.append(tmp)
        all_dic["|".join(i)]=np.std(tmp)
        all_std.append(np.std(tmp))
    min_std=min(all_std)
    #print(all_dic,min_std)
    for i in all_dic.keys():
        if all_dic[i]==min_std:
            #print("choose_min_std",i.split("|"))
            return i.split("|")
        
def choose_max_middle_loop(select_tmp):
    if len(select_tmp)==1:
        return select_tmp
    max_mid_stem=0
    all_dic={}
    for i in select_tmp:
        len_tmp=int(i[2].split("-")[0])-int(i[1].split("-")[1])
        if max_mid_stem<len_tmp:
            max_mid_stem=len_tmp
        if len_tmp not in all_dic.keys():
            all_dic[len_tmp]=[]
        all_dic[len_tmp].append(i)
    #for i in all_dic.keys():
        #print("choose_max_middle_loop",i, all_dic[i])
    #print(max_mid_stem,all_dic[max_mid_stem])
    return all_dic[max_mid_stem]

#checked
def check_suitable_stem_length(iM,stem_short,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    Cmotif = "[cC]{"+str(stem_short)+"}\w{"+str(loop1_short)+","+str(loop1_long)+"}[cC]{"+str(stem_short)+"}\w{"+str(loop2_short)+","+str(loop2_long)+"}[cC]{"+str(stem_short)+"}\w{"+str(loop3_short)+","+str(loop3_long)+"}[cC]{"+str(stem_short)+"}"
    #print(Cmotif)
    lis=[]
    resC = regex.finditer(Cmotif, iM)
    for i in resC:
        lis.append(i)
    if len(lis)>0:
        return True
    else:
        return False

def filter_all_comformation(confor_lis,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    fil=[]
    for i in confor_lis:
        stem1=i[0].split("-")
        stem2=i[1].split("-")
        stem3=i[2].split("-")
        stem4=i[3].split("-")
        loop1=int(stem2[0])-int(stem1[1])
        loop2=int(stem3[0])-int(stem2[1])
        loop3=int(stem4[0])-int(stem3[1])
        if loop1>=loop1_short and loop1<=loop1_long:
            if loop2>=loop2_short and loop2<=loop2_long:
                if loop3>=loop3_short and loop3<=loop3_long:
                    fil.append(i)
    return fil
    
def digest_conformation(iM,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    Ctract = "[cC]{"+str(stem_short)+",}"
    resC = regex.finditer(Ctract, iM)
    tract=[]
    for i in resC:
        s, e=i.span()[0], i.span()[1]
        tract.append(e-s)
    max_tract=max(tract)
    max_stem=min(max_tract,stem_long)
    all_conformation=[]
    loop_short=min(loop1_short,loop2_short,loop3_short)
    loop_long=max(loop1_long,loop2_long,loop3_long)
    for i in range(stem_short,max_stem+1):
        if check_suitable_stem_length(iM,i,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
            conformations=identifier_imotif_directedGraphic_alliM(iM, i, "Cc",loop_short,loop_long)
            all_conformation=all_conformation+conformations
            
    #print(len(all_conformation))
    return all_conformation

def digest_iM_pattern_regularExpre(stem_short,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    loop_short=min(loop1_short,loop2_short,loop3_short)
    loop_long=max(loop1_long,loop2_long,loop3_long)
    Cmotif = "[cC]{"+str(stem_short)+",}\w{"+str(loop_short)+","+str(loop_long)+"}[cC]{"+str(stem_short)+",}\w{"+str(loop_short)+","+str(loop_long)+"}[cC]{"+str(stem_short)+",}\w{"+str(loop_short)+","+str(loop_long)+"}[cC]{"+str(stem_short)+",}"
    return Cmotif

def combine_overlapped_iM_segment(seq,Cmotif):
    res = regex.finditer(Cmotif, seq, overlapped=True)
    iM_tmp_C=[]
    for kk in res:
        s, e = kk.span()[0], kk.span()[1]
        iM_tmp_C.append([s,e])
    if len(iM_tmp_C)==0:
        return iM_tmp_C
    iM_segment_C,non_iM_segment_C=combined_iMregion_and_NoniMregion(iM_tmp_C,[0,len(seq)])
    return iM_segment_C

def sorted_all_conformation_output_longest_shortest_each_position(greedy,conformation):
    dic={}
    for i in conformation:
        start=i[0].split("-")[0]
        if start not in dic.keys():
            dic[start]=[]
        dic[start].append(i)
    result={}
    if greedy==True:
        for i in dic.keys():
            result[i]=choose_longest_conformation(dic[i])
    else:
        for i in dic.keys():
            result[i]=choose_shortest_conformation(dic[i])
    return result
    
def choose_nonOverlapped_segment(sort):
    choosed=[str(sort[0][0])]
    sig=sort[0][0]+sort[0][1]
    end=sort[-1][0]+sort[-1][1]
    while sig <= end:
        if sig > sort[-1][0]:
            break
        for i in sort:
            if i[0]>=sig:
                sig=i[0]+i[1]
                choosed.append(str(i[0]))
                continue
    return choosed

def digest_nonOverlapped_segment(conformation_each_position):
    sort=[]
    for i in conformation_each_position.keys():
        rep=conformation_each_position[i][0]
        length=int(rep[-1].split("-")[-1])-int(rep[0].split("-")[0])
        sort.append([int(i),length])
    sort=merge_sort(sort,0)
    choosed=choose_nonOverlapped_segment(sort)
    new_conformation_each_position={}
    for i in conformation_each_position.keys():
        if i in choosed:
            new_conformation_each_position[i]=conformation_each_position[i]
    average_set,mid_loop_long_set=digest_Overlapped_segment(new_conformation_each_position)
    return average_set,mid_loop_long_set
            
def digest_Overlapped_segment(conformation_each_position):
    average_set=[]
    mid_loop_long_set=[]
    
    for i in conformation_each_position.keys():
        max_stem_set_tmp=conformation_each_position[i]
        
        average_set_tmp=choose_min_std(max_stem_set_tmp)
        average_set.append(average_set_tmp)
        
        mid_loop_long_set_tmp=choose_max_middle_loop(max_stem_set_tmp)
        mid_loop_long_set_long_mid_tmp=choose_min_std(mid_loop_long_set_tmp)
        mid_loop_long_set.append(mid_loop_long_set_long_mid_tmp)

    return average_set,mid_loop_long_set

def check_nuc_requirement(seq,conformation):
    new=[]
    for i in conformation:
        s=int(i[0].split("-")[0])
        e=int(i[-1].split("-")[-1])
        if check_iM(seq[s:e]):
            new.append(i)
    return new

def digest_iM(seq, overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    C_motif= digest_iM_pattern_regularExpre(stem_short,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)

    iM_segment=combine_overlapped_iM_segment(seq,C_motif)

    if len(iM_segment)==0:
        return iM_segment #没有潜在的iM

    all_result=[]
   
    for i in iM_segment:
        start,end=i[0],i[1]
        iM_combine_seg=seq[start:end]
        all_conformation=digest_conformation(iM_combine_seg,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
        if all_conformation==[]:
            #print("no iM in the segment",seq)
            continue
        else:
            new_conformation_fil=filter_all_comformation(all_conformation,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
            new_conformation=check_nuc_requirement(iM_combine_seg,new_conformation_fil)
            if new_conformation==[]:
                #print("no iM in the segment fitting requirement",seq)
                continue
            else:
                conformation_each_position=sorted_all_conformation_output_longest_shortest_each_position(greedy,new_conformation)
                if overlapped==True:
                    average_set,mid_loop_long_set=digest_Overlapped_segment(conformation_each_position)
                else:
                    average_set,mid_loop_long_set=digest_nonOverlapped_segment(conformation_each_position)
                if len(average_set) != len(mid_loop_long_set):
                    print("Wrong in digest_iM")
                for i in range(len(average_set)):
                    all_result.append([start,end,iM_combine_seg,average_set[i],mid_loop_long_set[i]])
    return all_result

def digest_iM_with_all_conformation(seq, overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    C_motif= digest_iM_pattern_regularExpre(stem_short,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)

    iM_segment=combine_overlapped_iM_segment(seq,C_motif)

    if len(iM_segment)==0:
        return iM_segment,[] #没有潜在的iM

    all_result=[]
    all_conformation_result=[]
    for i in iM_segment:
        start,end=i[0],i[1]
        iM_combine_seg=seq[start:end]
        all_conformation=digest_conformation(iM_combine_seg,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
        if all_conformation==[]:
            #print("no iM in the segment",seq)
            continue
        else:
            new_conformation_fil=filter_all_comformation(all_conformation,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
            new_conformation=check_nuc_requirement(iM_combine_seg,new_conformation_fil)
            #print(len(new_conformation))
            if new_conformation==[]:
                #print("no iM in the segment fitting requirement",seq)
                continue
            else:
                conformation_each_position=sorted_all_conformation_output_longest_shortest_each_position(greedy,new_conformation)
                if overlapped==True:
                    average_set,mid_loop_long_set=digest_Overlapped_segment(conformation_each_position)
                else:
                    average_set,mid_loop_long_set=digest_nonOverlapped_segment(conformation_each_position)
                if len(average_set) != len(mid_loop_long_set):
                    print("Wrong in digest_iM")
                for i in range(len(average_set)):
                    all_result.append([start,end,iM_combine_seg,average_set[i],mid_loop_long_set[i]])
                all_conformation_result.append([start,end,iM_combine_seg,new_conformation])
    
    return all_result,all_conformation_result



def detect_excat_7_segment(iM,conformation):
    stem1=[int(i) for i in conformation[0].split("-")]
    stem2=[int(i) for i in conformation[1].split("-")]
    stem3=[int(i) for i in conformation[2].split("-")]
    stem4=[int(i) for i in conformation[3].split("-")]
    stem_length=stem1[1]-stem1[0]
    loop1_length=int(stem2[0])-int(stem1[1])
    loop2_length=int(stem3[0])-int(stem2[1])
    loop3_length=int(stem4[0])-int(stem3[1])
    if stem_length*4+loop1_length+loop2_length+loop3_length!=stem4[1]-stem1[0]:
        print("False in detect_excat_region")
    stem1_seg=iM[stem1[0]:stem1[1]]
    loop1_seg=iM[stem1[1]:stem2[0]]
    stem2_seg=iM[stem2[0]:stem2[1]]
    loop2_seg=iM[stem2[1]:stem3[0]]
    stem3_seg=iM[stem3[0]:stem3[1]]
    loop3_seg=iM[stem3[1]:stem4[0]]
    stem4_seg=iM[stem4[0]:stem4[1]]
    if iM[stem1[0]:stem4[1]] != stem1_seg+loop1_seg+stem2_seg+loop2_seg+stem3_seg+loop3_seg+stem4_seg:
        print("Wrong in detect_excat_7_segment")
    return stem1[0], stem4[1],[stem1_seg,loop1_seg,stem2_seg,loop2_seg,stem3_seg,loop3_seg,stem4_seg]

def iM_list_operation(seq,iM_list,sig):
    result=[]
    for i in iM_list:
        pim_seq=i[2]
        pim_con=i[sig]
        start_tmp, end_tmp, segment_lis=detect_excat_7_segment(pim_seq,pim_con)
        start_true, end_true=i[0]+start_tmp, i[0]+end_tmp
        start_show, end_show=i[0]+start_tmp+1, i[0]+end_tmp
        if seq[start_true:end_true] != "".join(segment_lis):
            print(seq[start_true:end_true],"".join(segment_lis))
            print("wrong in digest_iM_for_conformation_choosed")
        tmp_lis=[start_show, end_show]+segment_lis
        result.append(tmp_lis)
    return result

def all_conformation_list_operation(seq,all_conformation_list):
    result=[]
    for i in all_conformation_list:
        pim_seq=i[2]
        for j in i[-1]:
            pim_con=j
            start_tmp, end_tmp, segment_lis=detect_excat_7_segment(pim_seq,j)
            start_true, end_true=i[0]+start_tmp, i[0]+end_tmp
            start_show, end_show=i[0]+start_tmp+1, i[0]+end_tmp
            if seq[start_true:end_true] != "".join(segment_lis):
                print(seq[start_true:end_true],"".join(segment_lis))
                print("wrong in digest_iM_for_conformation_choosed")
            tmp_lis=[start_show, end_show]+segment_lis
            result.append(tmp_lis)
    return result

def digest_iM_for_webServer(seq, conformation_choosed,overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    iM_list=digest_iM(seq, overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
    if len(iM_list)==0:
        return []
    if conformation_choosed==1:
        sig=3
    elif conformation_choosed==2:
        sig=4
    else:
        print("Wrong conformation choosed")
        return []
    result=iM_list_operation(seq,iM_list,sig)
    return result

def digest_iM_for_software(seq, all_conformation_out,conformation_choosed,overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long):
    iM_list,all_conformation_list=digest_iM_with_all_conformation(seq, overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
    if len(iM_list)==0:
        return False
    dic={}
    if conformation_choosed==1:
        result=iM_list_operation(seq,iM_list,3)
        dic["form1"]=result
    elif conformation_choosed==2:
        result=iM_list_operation(seq,iM_list,4)
        dic["form2"]=result
    elif conformation_choosed==3:
        result1=iM_list_operation(seq,iM_list,3)
        dic["form1"]=result1
        result2=iM_list_operation(seq,iM_list,4)
        dic["form2"]=result2
    else:
        print("Wrong conformation choosed")
        return False
    
    if all_conformation_out==True:
        conformation_result=all_conformation_list_operation(seq,all_conformation_list)
        dic["all"]=conformation_result
    return dic



def get_reversed_chain(seq):
    dic={
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C"}
    new_seq=""
    for i in seq:
        if i not in dic.keys():
            new_seq="N"+new_seq
        else:
            new_seq=dic[i]+new_seq
    return new_seq
def get_complementary_chain(seq):
    dic={
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C"}
    new_seq=""
    for i in seq:
        if i not in dic.keys():
            new_seq=new_seq+"N"
        else:
            new_seq=new_seq+dic[i]
    return new_seq



def quantify_GQS(seq):
  Ctract = "[cC]{1,}"
  Gtract = "[gG]{1,}"
  score=0
  resC = regex.finditer(Ctract, seq)
  resG = regex.finditer(Gtract, seq)
  for i in resG:
    s, e=i.span()[0], i.span()[1]
    if e-s<4:
        score=score+(e-s)**2
    else:
        score=score+(e-s)*4
  for i in resC:
    s, e=i.span()[0], i.span()[1]
    if e-s<4:
        score=score-(e-s)**2
    else:
        score=score-(e-s)*4
  #score=max([score*stemlength/looplength,0])
  #score=score*stemlength/looplength
  #score=max([score*stemlength/looplength/4,0])
  #score=score*stemlength/looplength/4
  return score

def digest_occur(seq,pattern,length):
    resC = regex.finditer(pattern, seq, overlapped=True)
    tract=[]
    for i in resC:
        s, e=i.span()[0], i.span()[1]
        tract.append([s,e])
    return len(tract)/length


def select_loop(iM_lis):
    tmp = []
    for j in range(len(iM_lis)-1):
        tmp.append(int(iM_lis[j+1].split("-")[0])-int(iM_lis[j].split("-")[1]))
    return max(tmp), min(tmp), np.mean(tmp), np.std(tmp)


def getFeatures(ctrac1,loop1,ctrac2,loop2,ctrac3,loop3,ctrac4):
    ctracNum=len(ctrac1)
    imotif_seq=ctrac1+loop1+ctrac2+loop2+ctrac3+loop3+ctrac4
    imotif_len=len(ctrac1)+len(ctrac2)+len(ctrac3)+len(ctrac4)+len(loop1)+len(loop2)+len(loop3)

    all_A_density=imotif_seq.count("A")/imotif_len
    all_C_density=imotif_seq.count("C")/imotif_len
    all_G_density=imotif_seq.count("G")/imotif_len
    all_T_density=imotif_seq.count("T")/imotif_len

    allLoop_seq=loop1+loop2+loop3
    allLoop_len=len(allLoop_seq)
    allLoop_A_density=allLoop_seq.count("A")/allLoop_len
    allLoop_C_density=allLoop_seq.count("C")/allLoop_len
    allLoop_G_density=allLoop_seq.count("G")/allLoop_len
    allLoop_T_density=allLoop_seq.count("T")/allLoop_len

    loop2_A_density=loop2.count("A")/len(loop2)
    loop2_C_density=loop2.count("C")/len(loop2)
    loop2_G_density=loop2.count("G")/len(loop2)
    loop2_T_density=loop2.count("T")/len(loop2)
    
    side_loop_A_density=(loop1.count("A")+loop3.count("A"))/(len(loop1)+len(loop3))
    side_loop_C_density=(loop1.count("C")+loop3.count("C"))/(len(loop1)+len(loop3))
    side_loop_G_density=(loop1.count("G")+loop3.count("G"))/(len(loop1)+len(loop3))
    side_loop_T_density=(loop1.count("T")+loop3.count("T"))/(len(loop1)+len(loop3))

    loop1_len=len(loop1)
    loop2_len=len(loop2)
    loop3_len=len(loop3)
    maxLoop=max(len(loop1),len(loop2),len(loop3))
    minLoop=min(len(loop1),len(loop2),len(loop3))
    
    if loop1_len>loop3_len:
        max_side_loop=loop1
        min_side_loop=loop3
    else:
        max_side_loop=loop3
        min_side_loop=loop1
        
    max_side_loop_A_density=max_side_loop.count("A")/len(max_side_loop)
    max_side_loop_C_density=max_side_loop.count("C")/len(max_side_loop)
    max_side_loop_G_density=max_side_loop.count("G")/len(max_side_loop)
    max_side_loop_T_density=max_side_loop.count("T")/len(max_side_loop)
    min_side_loop_A_density=min_side_loop.count("A")/len(min_side_loop)
    min_side_loop_C_density=min_side_loop.count("C")/len(min_side_loop)
    min_side_loop_G_density=min_side_loop.count("G")/len(min_side_loop)
    min_side_loop_T_density=min_side_loop.count("T")/len(min_side_loop)
    return [
        ctracNum,
        imotif_len,
        all_A_density,
        all_C_density,
        all_G_density,
        all_T_density,
        allLoop_len,
        loop2_len,
        loop1_len+loop3_len,
        max(loop1_len,loop3_len),
        min(loop1_len,loop3_len),
        maxLoop,
        minLoop,
        allLoop_A_density,
        allLoop_C_density,
        allLoop_G_density,
        allLoop_T_density,
        loop2_A_density,
        loop2_C_density,
        loop2_G_density,
        loop2_T_density,
        side_loop_A_density,
        side_loop_C_density,
        side_loop_G_density,
        side_loop_T_density,
        max_side_loop_A_density,
        max_side_loop_C_density,
        max_side_loop_G_density,
        max_side_loop_T_density,
        min_side_loop_A_density,
        min_side_loop_C_density,
        min_side_loop_G_density,
        min_side_loop_T_density
    ]


def read_iM_putative_result(file):
    dic={}
    with open(file,"r") as f:
        for line in f:
            seq=line.strip().split("\t")
            if seq[0] !="ID":
                name="|".join(seq)
                dic[name]=[seq[7],seq[8],seq[9],seq[10],seq[11],seq[12],seq[13]]
    return dic
                
            
def main():
    parser = argparse.ArgumentParser(description="Welcome to iM-Seeker.This is to predict DNA i-motif folding status and folding strength")
    parser.add_argument("--sequence", help="Show the path of sequence file in fasta format, default=None", type=str, required=True)
    parser.add_argument("--classification_model", help="Show the path of classification model , default=None", type=str, required=True)
    parser.add_argument("--regression_model", help="Show the path of regression model, default=None", type=str, required=True)
    parser.add_argument("--overlapped", help="Choose overlapped strategy (type 1) or non-overlapped strategy (type 2), default=2", type=int)
    parser.add_argument("--greedy", help="Choose greedy strategy (type 1) or non-greedy strategy (type 2), default=1", type=int)
    parser.add_argument("--stem_short", help="Set the lower boundary of C-tract length (can not be lower than 3), default=3", type=int)
    parser.add_argument("--stem_long", help="Set the higher boundary of C-tract length, default=100", type=int)
    parser.add_argument("--loop1_short", help="Set the lower boundary of loop1 length (can not be lower than 1, we recommand at least 1), default=1", type=int)
    parser.add_argument("--loop1_long", help="Set the higher boundary of loop1 length, default=12", type=int)
    parser.add_argument("--loop2_short", help="Set the lower boundary of loop2 length (can not be lower than 1,we recommand at least 1), default=1", type=int)
    parser.add_argument("--loop2_long", help="Set the higher boundary of loop2 length, default=12", type=int)
    parser.add_argument("--loop3_short", help="Set the lower boundary of loop3 length (can not be lower than 1,we recommand at least 1), default=1", type=int)
    parser.add_argument("--loop3_long", help="Set the higher boundary of loop3 length, default=12", type=int)
    parser.add_argument("--representative_conformation", help="Set the output of representative conformation (1 for average, 2 for shorter side loop), default=1", type=int)
    parser.add_argument("--output_folder", help="Show the path of output folder , default=folder containing input fasta file", type=str)

    parser.set_defaults(
        overlapped=2,
        greedy=1,
        stem_short=3,
        stem_long=100,
        loop1_short=1,
        loop1_long=12,
        loop2_short=1,
        loop2_long=12,     
        loop3_short=1,
        loop3_long=12,
        representative_conformation=1,
        output_folder=""
        )
    
    args = parser.parse_args()
    print(args)
    sequence=args.sequence
    try:
        seq=input_chr_file(sequence)
    except:
        print("Check the fasta file.")
        exit()

    overlapped=args.overlapped
    try:
        overlapped=int(overlapped)
    except:
        print("Check the overlapped parameter.")
        exit()
    if overlapped==1:
        overlapped=True
    elif overlapped==2:
        overlapped=False
    else:
        print("Check the overlapped parameter.")
        exit()
        
    greedy=args.greedy
    try:
        greedy=int(greedy)
    except:
        print("Check the greedy parameter.")
        exit()
    if greedy==1:
        greedy=True
    elif greedy==2:
        greedy=False
    else:
        print("Check the greedy parameter.")
        exit()
        
    stem_short=args.stem_short
    try:
        stem_short=int(stem_short)
    except:
        print("Check the stem_short parameter.")
        exit()
    if stem_short<3:
        print("Check the stem_short parameter.")
        exit()
        
    stem_long=args.stem_long
    try:
        stem_long=int(stem_long)
    except:
        print("Check the stem_long parameter.")
        exit()
    if stem_long<3 or stem_long<stem_short:
        print("Check the stem_long parameter.")
        exit()
        
    loop1_short=args.loop1_short
    try:
        loop1_short=int(loop1_short)
    except:
        print("Check the loop1_short parameter.")
        exit()
    if loop1_short<1:
        print("Check the loop1_short parameter.")
        exit()
        
    loop1_long=args.loop1_long
    try:
        loop1_long=int(loop1_long)
    except:
        print("Check the loop1_long parameter.")
        exit()
    if loop1_long<1 or loop1_long<loop1_short:
        print("Check the loop1_long parameter.")
        exit()

    loop2_short=args.loop2_short
    try:
        loop2_short=int(loop2_short)
    except:
        print("Check the loop2_short parameter.")
        exit()
    if loop2_short<1:
        print("Check the loop2_short parameter.")
        exit()
        
    loop2_long=args.loop2_long
    try:
        loop2_long=int(loop2_long)
    except:
        print("Check the loop2_long parameter.")
        exit()
    if loop2_long<1 or loop2_long<loop2_short:
        print("Check the loop2_long parameter.")
        exit()

    loop3_short=args.loop3_short
    try:
        loop3_short=int(loop3_short)
    except:
        print("Check the loop3_short parameter.")
        exit()
    if loop3_short<1:
        print("Check the loop3_short parameter.")
        exit()
        
    loop3_long=args.loop3_long
    try:
        loop3_long=int(loop3_long)
    except:
        print("Check the loop3_long parameter.")
        exit()
    if loop3_long<1 or loop3_long<loop3_short:
        print("Check the loop3_long parameter.")
        exit()
    
    representative_conformation=args.representative_conformation
    if representative_conformation not in [1,2]:
        print("Check the representative_conformation parameter.")
        exit()
    
    output_folder=args.output_folder
    if output_folder!="":
        if not os.path.isdir(output_folder):
            os.makedirs(output_folder, exist_ok=True)
            
    if output_folder=="":
        folder_path, file_name = os.path.split(sequence)
        output_folder=folder_path

    
    print("sequence:",sequence)
    print("classification_model",args.classification_model)
    print("regression_model",args.regression_model)
    print("overlapped:",overlapped)
    print("greedy:",greedy)
    print("stem_short:",stem_short)
    print("stem_long:",stem_long)
    print("loop1_short:",loop1_short)
    print("loop1_long:",loop1_long)
    print("loop2_short:",loop2_short)
    print("loop2_long:",loop2_long)
    print("loop3_short:",loop3_short)
    print("loop3_long:",loop3_long)
    print("representative_conformation:",representative_conformation)
    print("output_folder:",output_folder)

        
    rep1=os.path.join(output_folder,"iM-seeker_result_average_conformation.txt")
    rep2=os.path.join(output_folder,"iM-seeker_result_side_shorter_conformation.txt")
    
    if representative_conformation==1:
        out1=open(rep1,"w")
        out1.write("\t".join(["ID","start","end","stem_length","loop1_length","loop2_length","loop3_length"])+"\t")
        out1.write("\t".join(["stem1","loop1","stem2","loop2","stem3","loop3","stem4"])+"\n")
    if representative_conformation==2:
        out2=open(rep2,"w")
        out2.write("\t".join(["ID","start","end","stem_length","loop1_length","loop2_length","loop3_length"])+"\t")
        out2.write("\t".join(["stem1","loop1","stem2","loop2","stem3","loop3","stem4"])+"\n")
        

    for i in seq.keys():
        dic=digest_iM_for_software(seq[i], False,representative_conformation,overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
        if dic != False:
            if "form1" in dic.keys():
                for j in dic["form1"]:
                    out1.write("\t".join([i+"+",str(j[0]),str(j[1]),str(len(j[2])),str(len(j[3])),str(len(j[5])),str(len(j[7]))])+"\t"+"\t".join(j[2:])+"\n")
            if "form2" in dic.keys():
                for j in dic["form2"]:
                    out2.write("\t".join([i+"+",str(j[0]),str(j[1]),str(len(j[2])),str(len(j[3])),str(len(j[5])),str(len(j[7]))])+"\t"+"\t".join(j[2:])+"\n")
    for i in seq.keys():
        minus_seq=get_complementary_chain(seq[i])
        dic=digest_iM_for_software(minus_seq, False,representative_conformation,overlapped, greedy,stem_short,stem_long,loop1_short,loop1_long,loop2_short,loop2_long,loop3_short,loop3_long)
        if dic != False:
            if "form1" in dic.keys():
                for j in dic["form1"]:
                    out1.write("\t".join([i+"-",str(j[0]),str(j[1]),str(len(j[2])),str(len(j[3])),str(len(j[5])),str(len(j[7]))])+"\t"+"\t".join(j[2:])+"\n")
            if "form2" in dic.keys():
                for j in dic["form2"]:
                    out2.write("\t".join([i+"-",str(j[0]),str(j[1]),str(len(j[2])),str(len(j[3])),str(len(j[5])),str(len(j[7]))])+"\t"+"\t".join(j[2:])+"\n")
       
    if representative_conformation==1:
        out1.close()
        result_putative_iM=read_iM_putative_result(rep1)
    if representative_conformation==2:
        out2.close()
        result_putative_iM=read_iM_putative_result(rep2)

    classification_model = pickle.load(open(args.classification_model, 'rb'))
    regression_model = pickle.load(open(args.regression_model, 'rb'))
    
    out_file=os.path.join(output_folder,"iM-seeker_final_prediction.txt")
    out=open(out_file,"w")
    out.write("\t".join(["Putative i-motif id","Folding status","Folding strength","Folding probability"])+"\n")
    for i in result_putative_iM.keys():
        feature_list=getFeatures(result_putative_iM[i][0],result_putative_iM[i][1],result_putative_iM[i][2],result_putative_iM[i][3],
                                 result_putative_iM[i][4],result_putative_iM[i][5],result_putative_iM[i][6])
        classification_result=classification_model.predict([feature_list])
        regression_result=regression_model.predict([feature_list])
        out.write("\t".join([i,str(classification_result[0]),str(regression_result[0]),"%s" % list(classification_model.predict_proba([feature_list])[0])])+"\n")
    out.close()
if __name__ == "__main__":
    main()

