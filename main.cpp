#include <iostream>

#include "centrality.hpp"

int main(int argc, char** argv){

    if( argc != 4 && argc != 5){
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. MA2013\n");
        printf("<arg3> task, 1: MCCMN; 2: SCCMN; 3: PSCCMN. e.g. 3\n");
        printf("<arg4> thread number(optional), e.g. 10\n");
        exit(0);
    }
    cout << "Hello, World!" << endl;
    string source_path = "/Users/zhouxj/Documents/1-Research/Datasets";
    Graph g;
    g.threadnum = 10;
    g.dataset = "TestMN";//MA2013
    int task = 4;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1]: " << argv[1] << endl;//source path
        source_path = argv[1];
        if(argc > 2){
            cout << "argv[2]: " << argv[2] << endl;//dataset
            g.dataset = argv[2];
        }
        if(argc > 3){
            cout << "argv[3]: " << argv[3] << endl;//dataset
            task = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4]: " << argv[4] << endl;//dataset
            g.threadnum = stoi(argv[4]);
        }
    }

    cout<<"Dataset: "<<g.dataset<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;

    g.graph_path = source_path + "/" + g.dataset + "/" + g.dataset;
    g.ReadMultiGraph(g.graph_path);
    g.Preprocess();

    switch (task) {
        case 1:{
            g.CentralityCompute1();// Baseline
            break;
        }
        case 2:{
            g.CentralityCompute2();// Algorithm 1
            break;
        }
        case 3:{
            g.CentralityCompute3();// Algorithm 2
            break;
        }
        case 4:{
            g.TopKCentrality(10);
            break;
        }
        case 5:{
            g.ActivityProcess();
            break;
        }
        default:
            break;
    }


    return 0;
}
