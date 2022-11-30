//
// Created by Xinjie ZHOU on 27/11/2022.
//
#ifndef CENTRALITY_HPP
#define CENTRALITY_HPP

#include "centrality.h"

void Graph::ReadMultiGraph(string filename){
    string r_edges = filename + "_multiplex.edges";

    ifstream inFile(r_edges, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << r_edges << endl;
        exit(1);
    }

    string line;
    int layerID;
    int ID1,ID2,weight;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> layer_num >> node_num)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
    cout<<"Layer number: "<<layer_num<<" ; Node number: "<<node_num<<endl;
    multiGraph.assign(layer_num,vector<unordered_map<int,int>>(node_num));
    degrees.assign(layer_num,vector<int>(node_num));
    nodeSets.assign(layer_num,unordered_set<int>());
    layerNumbers.assign(layer_num, make_pair(0,0));

    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        if (!(iss >> layerID >> ID1 >> ID2 >> weight)){
            cout<<"Wrong input syntax!"<<endl;
            exit(1);
        }

        layerID-=1; ID1-=1; ID2-=1;
        multiGraph[layerID][ID1].insert({ID2,weight});
        multiGraph[layerID][ID2].insert({ID1,weight});

    }
    inFile.close();
    for(int id=0;id<node_num;++id){
        for(layerID=0;layerID<layer_num;++layerID){
            degrees[layerID][id] = multiGraph[layerID][id].size();
            if(!multiGraph[layerID][id].empty()){//if not empty
                layerNumbers[layerID].first++;
                layerNumbers[layerID].second+=degrees[layerID][id];
                nodeSets[layerID].insert(id);
            }
        }
    }
    for(layerID=0;layerID<layer_num;++layerID){
        cout<<"Layer "<<layerID<<": "<<layerNumbers[layerID].first<<" "<<layerNumbers[layerID].second<<endl;
    }
}

void Graph::Preprocess() {

    nodeSetsLCC.assign(layer_num,set<int>());
    vector<vector<int>> CCs;
    for(int layerID=0;layerID<layer_num;++layerID){
        cout<<"Layer "<<layerID<<". ";
        CCs.emplace_back(DFS_CC(multiGraph[layerID],nodeSets[layerID],nodeSetsLCC[layerID],node_num));
    }
    //write to disk
    string filename = graph_path + "_CCs.txt";
    ofstream outFile(filename);
    //line 1:
    outFile<<CCs.size()<<endl;
    for(int i=0;i<CCs.size();++i){
        outFile<<CCs[i].size()<<endl;
        for(int j=0;j<CCs[i].size();++j){
            outFile<<CCs[i][j]<<" ";
        }
        outFile<<endl;
    }
    outFile.close();
}

//function of checking the connectivity
vector<int> Graph::DFS_CC(vector<unordered_map<int,int>> & Edges, unordered_set<int> & set_A, set<int> & set_LCC, int nodenum) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    unordered_set<int> set_B;//nodes visited for current component
    set_B.clear();
    int item_id,temp_id;
    vector<bool> flag_visited(nodenum,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<unordered_set<int>,int> LCC;
    vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<nodenum;++i){
//        set_A.insert(i);
//    }
    int seed = *set_A.begin();
    stack_A.push(seed);
    set_A.erase(seed);
    set_B.insert(seed);
    flag_visited[seed] = true;
    //Get the connected components by DFS
    while(!set_A.empty()) {//if not finish
        temp_num = 0;
        while (!stack_A.empty()) {
            item_id = stack_A.top();
            stack_A.pop();
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                temp_num += 1;
                if (!flag_visited[temp_id]) {//if not visited
                    stack_A.push(temp_id);
                    set_A.erase(temp_id);
                    set_B.insert(temp_id);
                    flag_visited[temp_id] = true;
                }
            }
        }
        if (set_B.size() > LCC.first.size()) {
            LCC.first.clear();
            LCC.first = set_B;
            LCC.second = temp_num;// /2
        }
        assert(!set_B.empty());
        CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
        ++component_i;
        set_B.clear();
        if (!set_A.empty()) {
            stack_A.push(*set_A.begin());
            set_B.insert(*set_A.begin());
            flag_visited[*set_A.begin()] = true;
            set_A.erase(*set_A.begin());
        } else {
            break;
        }
    }
    if(component_i==1){
        cout<<"This graph has only one connected component. ";
        cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
        cout<<"Edges size of graph: "<< LCC.second << endl;
    }else{
        cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
    }
    for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
        set_LCC.insert(*it);
    }
    std::sort(CCs.begin(), CCs.end());
    return CCs;
//    return component_i;
}

//baseline algorithm
void Graph::CentralityCompute1() {
    cout<<"This is testing for CCMN."<<endl;

    betweennessScores.assign(layer_num,vector<double>(node_num,0.0));
    closenessScores.assign(layer_num,vector<double>(node_num,0.0));

    Timer tt;
    tt.start();

    for(int layerID=0;layerID<layer_num;++layerID){
        cout<<"Layer "<<layerID<<endl;
        Betweenness(multiGraph[layerID],nodeSetsLCC[layerID],betweennessScores[layerID]);
        Closeness(multiGraph[layerID],nodeSetsLCC[layerID],closenessScores[layerID]);
    }

    tt.stop();
    cout<<"Overall time for CCMN: "<<tt.GetRuntime()<<" s."<<endl;
    WriteResults(graph_path+"_CC1.txt");
}
//function of algorithm 1
void Graph::CentralityCompute2() {
    cout<<"This is testing for SCCMN."<<endl;

    betweennessScores.assign(layer_num,vector<double>(node_num,0.0));
    closenessScores.assign(layer_num,vector<double>(node_num,0.0));

    Timer tt;
    tt.start();

    for(int layerID=0;layerID<layer_num;++layerID){
        cout<<"Layer "<<layerID<<endl;
        Centrality(multiGraph[layerID],nodeSetsLCC[layerID],layerID);
    }

    tt.stop();
    cout<<"Overall time for SCCMN: "<<tt.GetRuntime()<<" s."<<endl;
    WriteResults(graph_path+"_CC2.txt");
}

//function of algorithm 2
void Graph::CentralityCompute3() {
    cout<<"This is testing for PSCCMN."<<endl;

    betweennessScores.assign(layer_num,vector<double>(node_num,0.0));
    closenessScores.assign(layer_num,vector<double>(node_num,0.0));

    Timer tt;
    tt.start();

    boost::thread_group threadf;
    for(int layerID=0;layerID<layer_num;++layerID){
        cout<<"Layer "<<layerID<<endl;
        threadf.add_thread(new boost::thread(&Graph::CentralityParallel, this, boost::ref(multiGraph[layerID]), boost::ref(nodeSetsLCC[layerID]), layerID));
    }
    threadf.join_all();

    tt.stop();
    cout<<"Overall time for PSCCMN: "<<tt.GetRuntime()<<" s."<<endl;
    WriteResults(graph_path+"_CC3.txt");
}

//function of computing degree, betweenness, and closeness, simultaneously
void Graph::CentralityParallel(vector<unordered_map<int, int>> &Edges, set<int> &Nodes, int layerID) {
    vector<int> vNodes;
    for(auto it=Nodes.begin();it!=Nodes.end();++it){
        vNodes.emplace_back(*it);
    }
    int nodeNumSet = Nodes.size();
    int totalNum = 0;
    for(int i=0;i<layer_num;++i){
        totalNum += nodeSetsLCC[i].size();
    }
//    int threadNum = threadnum/layer_num;//thread number for this layer
    int threadNum = 1;
    threadNum = max(threadNum,nodeNumSet*threadnum/totalNum);//thread number for this layer
    cout<<"Thread number for layer "<<layerID<<": "<<threadNum<<endl;

    if(nodeNumSet>threadNum){
        int step=nodeNumSet/threadNum;
        boost::thread_group threadf;
        for(int i=0;i<threadNum;i++){
            pair<int,int> p;
            p.first=i*step;
            if(i==threadNum-1)
                p.second=nodeNumSet;
            else
                p.second=(i+1)*step;
            threadf.add_thread(new boost::thread(&Graph::CentralityP, this, p, boost::ref(multiGraph[layerID]), boost::ref(vNodes), layerID));
        }
        threadf.join_all();
    }else{
        boost::thread_group threadf;
        for(int pid=0;pid<nodeNumSet;++pid) {
            threadf.add_thread(new boost::thread(&Graph::CentralityP, this, make_pair(pid,pid+1), boost::ref(multiGraph[layerID]), boost::ref(vNodes), layerID));
        }
        threadf.join_all();
    }

    cout<<"Done."<<endl;
}

void Graph::CentralityP(pair<int,int> pidRange, vector<unordered_map<int, int>> &Edges, vector<int> &Nodes, int layerID) {
    int nodeNum = Edges.size();
    int nodeNumSet = Nodes.size();

    for(int pid=pidRange.first;pid<pidRange.second;++pid) {
        //single-source shortest path for each vertex
        int node_start = Nodes[pid];

        benchmark::heap<2, int, int> pqueue(nodeNum);
        stack<int> S;
        vector<int> spNums(node_num, 0);//the number of the shortest paths from node_start to v
        vector<int> counts(node_num, 0);
        int item_id, temp_id;
        int item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<int> cost(node_num, INF);   //vector of cost
        vector<unordered_set<int>> pre(node_num);       //list of predecessor id
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.update(node_start, 0);
        counts[node_start] = 1;
        int count = 0;
        double CC = 0;
        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            pqueue.extract_min(item_id, item_dis);// top and delete min item
            S.push(item_id);
            spNums[item_id] += counts[item_id];
            CC += item_dis;
            count++;
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + it->second;
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pre[temp_id].clear();
                    pre[temp_id].insert(item_id);
                    pqueue.update(temp_id, temp_dis);
                    counts[temp_id] = spNums[item_id];
                } else if (cost[temp_id] == temp_dis) {
                    counts[temp_id] += spNums[item_id];
                    pre[temp_id].insert(item_id);
                }

            }
            closed[item_id] = true;
        }

        assert(count == nodeNumSet);
        CC = (nodeNumSet - 1) / CC;
        closenessScores[layerID][node_start] = CC;

        vector<double> dependency(node_num, 0.0);
        double temp = 0;
        while (!S.empty()) {
            item_id = S.top();
            S.pop();
            temp = (double) (1 + dependency[item_id]) / spNums[item_id];
            for (auto it = pre[item_id].begin(); it != pre[item_id].end(); ++it) {
                temp_id = *it;//predecessor
                dependency[temp_id] += temp * spNums[temp_id];
            }
            if (item_id != node_start) {
                betweennessScores[layerID][item_id] += dependency[item_id];
            }
        }
    }

}

//function of computing degree, betweenness, and closeness, simultaneously
void Graph::Centrality(vector<unordered_map<int, int>> &Edges, set<int> &Nodes, int layerID) {
    int nodeNum = Edges.size();
    int nodeNumSet = Nodes.size();

    //single-source shortest path for each vertex
    for(auto itd=Nodes.begin();itd!=Nodes.end();++itd){
        int node_start = *itd;

        benchmark::heap<2, int, int> pqueue(nodeNum);
        stack<int> S;
        vector<int> spNums(node_num,0);//the number of the shortest paths from node_start to v
        vector<int> counts(node_num,0);
        int item_id, temp_id;
        int item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<int> cost(node_num, INF);   //vector of cost
        vector<unordered_set<int>> pre(node_num);       //list of predecessor id
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.update(node_start, 0);
        counts[node_start]=1;
        int count = 0;
        double CC = 0;
        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            pqueue.extract_min(item_id, item_dis);// top and delete min item
            S.push(item_id);
            spNums[item_id] += counts[item_id];
            CC += item_dis;
            count++;
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + it->second;
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pre[temp_id].clear();
                    pre[temp_id].insert(item_id);
                    pqueue.update(temp_id, temp_dis);
                    counts[temp_id] = spNums[item_id];
                }else if(cost[temp_id] == temp_dis){
                    counts[temp_id] += spNums[item_id];
                    pre[temp_id].insert(item_id);
                }

            }
            closed[item_id] = true;
        }

        assert(count == nodeNumSet);
        CC = (nodeNumSet-1)/CC;
        closenessScores[layerID][node_start] = CC;

        vector<double> dependency(node_num,0.0);
        double temp = 0;
        while(!S.empty()){
            item_id = S.top();
            S.pop();
            temp = (double)(1+dependency[item_id])/spNums[item_id];
            for(auto it=pre[item_id].begin();it!=pre[item_id].end();++it){
                temp_id = *it;//predecessor
                dependency[temp_id] += temp * spNums[temp_id];
            }
            if(item_id != node_start){
                betweennessScores[layerID][item_id] += dependency[item_id];
            }
        }

    }
    cout<<"Done."<<endl;
}

//function of computing the shortest-path closeness centrality for all vertices
void Graph::Betweenness(vector<unordered_map<int, int>> &Edges, set<int> &Nodes, vector<double> &BCs) {
    int nodeNum = Edges.size();
    int nodeNumSet = Nodes.size();
//    BCs.assign(nodeNum,0.0);

    //single-source shortest path for each vertex
    for(auto itd=Nodes.begin();itd!=Nodes.end();++itd){
        int node_start = *itd;

        benchmark::heap<2, int, int> pqueue(nodeNum);
        stack<int> S;
        vector<int> spNums(node_num,0);//the number of the shortest paths from node_start to v
        vector<int> counts(node_num,0);
        int item_id, temp_id;
        int item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<int> cost(node_num, INF);   //vector of cost
        vector<unordered_set<int>> pre(node_num);       //list of predecessor id
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.update(node_start, 0);
        counts[node_start]=1;
        int count = 0;
        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            pqueue.extract_min(item_id, item_dis);// top and delete min item
            S.push(item_id);
            spNums[item_id] += counts[item_id];
            count++;
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + it->second;
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
                    pre[temp_id].clear();
                    pre[temp_id].insert(item_id);
                    pqueue.update(temp_id, temp_dis);
                    counts[temp_id] = spNums[item_id];
                }else if(cost[temp_id] == temp_dis){
                    counts[temp_id] += spNums[item_id];
                    pre[temp_id].insert(item_id);
                }

            }
            closed[item_id] = true;
        }

        assert(count == nodeNumSet);
        vector<double> dependency(node_num,0.0);
        double temp = 0;
        while(!S.empty()){
            item_id = S.top();
            S.pop();
            temp = (double)(1+dependency[item_id])/spNums[item_id];
            for(auto it=pre[item_id].begin();it!=pre[item_id].end();++it){
                temp_id = *it;//predecessor
                dependency[temp_id] += temp * spNums[temp_id];
            }
            if(item_id != node_start){
                BCs[item_id] += dependency[item_id];
            }
        }

    }
    cout<<"Done."<<endl;
}

//function of computing the shortest-path closeness centrality for all vertices
void Graph::Closeness(vector<unordered_map<int,int>> & Edges, set<int> & Nodes, vector<double> & CCs) {
    int nodeNum = Edges.size();
    int nodeNumSet = Nodes.size();
//    CCs.assign(nodeNum,0);

    //single-source shortest path for each vertex
    for(auto itd=Nodes.begin();itd!=Nodes.end();++itd){
        int node_start = *itd;

        benchmark::heap<2, int, int> pqueue(nodeNum);
        int item_id, temp_id;
        int item_dis, temp_dis;
        vector<bool> closed(node_num, false); //flag vector of whether closed
        vector<int> cost(node_num, INF);   //vector of cost
//        vector<int> pre(node_num, -1);       //vector of predecessor id
        double CC = 0;
        //Initiation of start node
        cost[node_start] = 0;//cost of start node
        pqueue.update(node_start, 0);
        int count = 0;
        //Iteration
        while (!pqueue.empty()) {//for every node in pqueue
            pqueue.extract_min(item_id, item_dis);// top and delete min item
            CC += item_dis;
            count++;
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                if (closed[temp_id])//if closed
                    continue;
                temp_dis = item_dis + it->second;
                if (cost[temp_id] > temp_dis) {//slack operation
                    cost[temp_id] = temp_dis;
//                    pre[temp_id] = item_id;
                    pqueue.update(temp_id, temp_dis);
                }
            }
            closed[item_id] = true;
        }
        assert(count == nodeNumSet);
        CC = (nodeNumSet-1)/CC;
        CCs[node_start] = CC;
    }
    cout<<"Done."<<endl;
}

void Graph::TopKCentrality(int k){
    ReadResults(graph_path+"_CC3.txt");
    vector<set<CompInt>> DegSets(layer_num,set<CompInt>());
    vector<set<CompDouble>> BCSets(layer_num,set<CompDouble>());
    vector<set<CompDouble>> CCSets(layer_num,set<CompDouble>());
    for(int layer=0;layer<layer_num;++layer){
        for(auto it=nodeSetsLCC[layer].begin();it!=nodeSetsLCC[layer].end();++it){
            DegSets[layer].insert(CompInt(*it,multiGraph[layer][*it].size()));
            BCSets[layer].insert(CompDouble(*it,betweennessScores[layer][*it]));
            CCSets[layer].insert(CompDouble(*it,closenessScores[layer][*it]));
        }
    }
    //top-10 highest degree vertex
    for(int layer=0;layer<layer_num;++layer){
        cout<<"Layer "<<layer<<endl;
        int id=0;
        vector<int> TopK;
        cout<<"ID";
        for(auto it=DegSets[layer].begin();id<k;++it){
            cout<<" "<<it->id;
            TopK.push_back(it->id); ++id;
        }
        cout<<"\nDeg";
        for(int i=0;i<TopK.size();++i){
            cout<<" "<<multiGraph[layer][TopK[i]].size();
        }
        cout<<"\nBC";
        for(int i=0;i<TopK.size();++i){
            cout<<" "<<betweennessScores[layer][TopK[i]];
        }
        cout<<"\nCC";
        for(int i=0;i<TopK.size();++i){
            cout<<" "<<closenessScores[layer][TopK[i]];
        }
        cout<<endl;
    }
    cout<<"Done."<<endl;
    //top-10 vertex
    vector<vector<vector<pair<int,double>>>> topK(layer_num,vector<vector<pair<int,double>>>(3,vector<pair<int,double>>()));
    for(int layer=0;layer<layer_num;++layer) {
        int count = 0;
        for(auto it=DegSets[layer].begin();count<k;++count,++it){
            topK[layer][0].emplace_back(it->id,it->value);
        }
        count = 0;
        for(auto it=BCSets[layer].begin();count<k;++count,++it){
            topK[layer][1].emplace_back(it->id,it->value);
        }
        count = 0;
        for(auto it=CCSets[layer].begin();count<k;++count,++it){
            topK[layer][2].emplace_back(it->id,it->value);
        }
    }

    for(int i=0;i<k;++i){
        for(int j=0;j<3;++j){
            for(int layer=0;layer<layer_num;++layer){
                cout<<topK[layer][j][i].first<<" ";
            }
        }
        cout<<endl;
    }
    cout<<"Done."<<endl;
}

void Graph::ActivityProcess(){
    ReadActivity(graph_path+"_activity.txt");

}

//function of writing the centrality result to disk
void Graph::WriteResults(string filename) {
    ofstream outFile(filename);
    //line 1:
    outFile<<layer_num<<" "<<node_num;
    for(int i=0;i<layer_num;++i){
        outFile<<" "<<nodeSetsLCC[i].size();
    }
    outFile<<endl;
    int ID;
    for(int layerID=0;layerID<layer_num;++layerID){
        for(auto it=nodeSetsLCC[layerID].begin();it!=nodeSetsLCC[layerID].end();++it){
            ID = *it;
            outFile<<layerID<<" "<<ID<<" "<<multiGraph[layerID][ID].size()<<" "<<betweennessScores[layerID][ID]<<" "<<closenessScores[layerID][ID]<<endl;
        }
    }

    outFile.close();
}

//function of reading the centrality result from disk
void Graph::ReadResults(string filename) {
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }

    string line;
    int layerNum1,layerNum2,layerNum3;
    int layerID;
    int ID,deg;

    getline(inFile,line);
    istringstream iss(line);
    if (!(iss >> layer_num >> node_num >> layerNum1 >> layerNum2 >> layerNum3)){
        cout<<"Wrong input syntax!"<<endl;
        exit(1);
    }
    cout<<"Layer number: "<<layer_num<<" ; Node number: "<<node_num<<" "<<layerNum1<<" "<<layerNum2<<" "<<layerNum3<<endl;

    betweennessScores.assign(layer_num,vector<double>(node_num,0.0));
    closenessScores.assign(layer_num,vector<double>(node_num,0.0));

    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        iss >> layerID >> ID >> deg;
        iss >> betweennessScores[layerID][ID] >> closenessScores[layerID][ID];

    }
    inFile.close();
}

void Graph::ReadActivity(string filename){
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }

    string line;

    int ID1,ID2;
    int time,time1,time2;
    string layerName;
    vector<vector<pair<int,int>>> activity;//<time, user number>
    vector<unordered_set<int>> users;
    activity.assign(layer_num,vector<pair<int,int>>());
    users.assign(layer_num,unordered_set<int>());

    while(getline(inFile,line)){
        if(line.empty())
            continue;
        istringstream iss(line);
        iss >> ID1 >> ID2 >> time >> layerName;
        if(layerName == "RT"){
            users[0].insert(ID1);
            users[0].insert(ID2);
            activity[0].emplace_back(time,users[0].size());
        }else if(layerName == "MT"){
            users[1].insert(ID1);
            users[1].insert(ID2);
            activity[1].emplace_back(time,users[1].size());
        }else if(layerName == "RE"){
            users[2].insert(ID1);
            users[2].insert(ID2);
            activity[2].emplace_back(time,users[2].size());
        }
    }
    inFile.close();

    time1 = INT32_MAX; time2 = 0;
    for(int i=0;i<layer_num;++i){
        time1 = min(time1,activity[i][0].first);
        time2 = max(time2,activity[i][activity[i].size()-1].first);
    }

    int hour = 3600;
    int spans = ceil((time2-time1)/hour);
    cout<<"Begin time: "<<time1<<" ; End time: "<<time2<<"; Time span: "<< spans <<" hours ("<<spans/24<<" days)."<<endl;
    //get the vector
    vector<vector<int>> volume(layer_num,vector<int>(spans,0));
    vector<vector<int>> userSum(layer_num,vector<int>(spans,0));
    int span = 0;

    for(int layerID=0;layerID<layer_num;++layerID){
        span = 0;
        for(int i=0;i<activity[layerID].size();++i){
            if(activity[layerID][i].first <= time1 + (span+1)*hour){
                volume[layerID][span]++;
                userSum[layerID][span] = activity[layerID][i].second;
            }else{
                span++;
                userSum[layerID][span] = userSum[layerID][span-1];
            }
        }
        assert(span == spans);
    }
    //print results
    cout<<"\nVolume"<<endl;
    for(int layerID=0;layerID<layer_num;++layerID){
        if(layerID == 0) cout<<"RT";
        else if(layerID == 1) cout<<"MT";
        else if(layerID == 2) cout<<"RE";
        for(int i=0;i<volume[layerID].size();++i){
            cout<<" "<<volume[layerID][i];
        }
        cout<<endl;
    }
    cout<<"\nUser"<<endl;
    for(int layerID=0;layerID<layer_num;++layerID){
        if(layerID == 0) cout<<"RT";
        else if(layerID == 1) cout<<"MT";
        else if(layerID == 2) cout<<"RE";
        for(int i=0;i<userSum[layerID].size();++i){
            cout<<" "<<userSum[layerID][i];
        }
        cout<<endl;
    }

    cout<<"Done."<<endl;
}

#endif //CENTRALITY_HPP