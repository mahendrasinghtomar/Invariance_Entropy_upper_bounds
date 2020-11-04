//
//  Helper_MST.h
//  SCOTSv02MST
//
//  Created by MST on 02.03.20.
//  Copyright © 2020 MST. All rights reserved.
//

#ifndef Helper_MST_h
#define Helper_MST_h

#include <sstream>
#include <fstream>
#include <string>
#include <unordered_map>

#include "scots.hh"

template<class F>
void print(const F &var, std::string name){
    /*To print array/vector.*/
	std::cout << name << "=(";
	for(int i=0;i<var.size();i++)
		std::cout << var[i] << ", ";
	std::cout << ")  ";
	return;
}

template< class F1, class F2, class F3, class state_type, class input_type, class WinDom>
void DeterminiseController_Xj_NPosts(int& int_ReachSet, F1& system_post, F2& radius_post, F3& reach_set, WinDom& win, scots::UniformGrid& ss, scots::UniformGrid& is, scots::Abstraction<state_type, input_type>& abs, const std::vector<scots::abs_type>& Bpartition, const scots::abs_type& Bsize, std::vector<std::vector<scots::abs_type>>& Xj, std::vector<scots::abs_type>& DControllerMinPost, std::vector<int>& B_NPosts) {
    // deterministic controller with such a input selected for which no. of posts is minimum
    std::vector<scots::abs_type> input_list;
    std::vector<state_type> post_list;
    std::vector<state_type> selInput_post_list;
    state_type xval;
    input_type uval;
    scots::abs_type selInput; // selected input
    int selInput_noPost; // number of post for the selected input
    for (int i = 0; i < Bsize; i++) {
        input_list = win.get_inputs(Bpartition[i]);
        for (int j = 0; j < input_list.size(); j++) {
            ss.itox(Bpartition[i], xval);
            is.itox(input_list[j], uval);
            if (int_ReachSet)
                post_list = abs.get_post(reach_set, xval, uval, "MST");
            else
                post_list = abs.get_post(system_post, radius_post, xval, uval);
            if (j == 0) {
                selInput_post_list = post_list;
                selInput_noPost = selInput_post_list.size();
                selInput = input_list[0];
            }
            else
            {
                if (post_list.size() < selInput_noPost) {
                    selInput = input_list[j];
                    selInput_post_list = post_list;
                    selInput_noPost = selInput_post_list.size();
                }
            }
        }
        DControllerMinPost[i] = selInput;
        B_NPosts[i] = selInput_noPost;
        for (int j = 0; j < selInput_noPost; j++)
            Xj[i].push_back(ss.xtoi(selInput_post_list[j]));
    }
}

template< class F1, class F2, class F3, class state_type, class input_type, class WinDom>
void minPost_controller(int& int_ReachSet, F1& system_post, F2& radius_post, F3& reach_set, WinDom& win, scots::UniformGrid& ss, scots::UniformGrid& is, scots::Abstraction<state_type, input_type>& abs, const std::vector<scots::abs_type>& Bpartition, const scots::abs_type& Bsize, std::vector<std::vector<scots::abs_type>>& controllerMinPost, std::vector<int>& B_NPosts) {
    // Controller with those inputs for which the number of posts is minimum.
    std::vector<scots::abs_type> input_list;
    std::vector<state_type> post_list;
    state_type xval;
    input_type uval;
    //scots::abs_type selInput; // selected input
    int min_noPost; // number of post for the selected input
    for (int i = 0; i < Bsize; i++) {
        input_list = win.get_inputs(Bpartition[i]);
        std::vector<scots::abs_type> minPost_input_list;
        for (int j = 0; j < input_list.size(); j++) {
            ss.itox(Bpartition[i], xval);
            is.itox(input_list[j], uval);
            if (int_ReachSet)
                post_list = abs.get_post(reach_set, xval, uval, "MST");
            else
                post_list = abs.get_post(system_post, radius_post, xval, uval);
            if (j == 0) {
                //selInput_post_list = post_list;
                min_noPost = post_list.size();
                //selInput = input_list[0];
                minPost_input_list.push_back(input_list[0]);
            }
            else
            {
                if (post_list.size() < min_noPost) {
                    minPost_input_list.clear(); 
                    //selInput = input_list[j];
                    //selInput_post_list = post_list;
                    min_noPost = post_list.size();
                    minPost_input_list.push_back(input_list[j]);
                }
                else if(post_list.size() == min_noPost)
                    minPost_input_list.push_back(input_list[j]);
            }
        }
        //DControllerMinPost[i] = selInput;
        controllerMinPost[i] = minPost_input_list;
        B_NPosts[i] = min_noPost;
        /*for (int j = 0; j < selInput_noPost; j++)
            Xj[i].push_back(ss.xtoi(selInput_post_list[j]));*/
    }
}


void compute_coarsePartition(const std::vector<scots::abs_type>& DControllerMinPost, 
                                    std::vector<int>& CoarsePartition, 
                                    std::vector<int>& Bpartition_coarsePno) {
    /*Coarse partition index taken to start from 1.
    CoarsePartition stores the unique control input values.
    Coarse partition index corresponding to a control input is = (the index of the control input in CoarsePartition) + 1. 
    Bpartition_coarsePno is expected to have been initialized to all zeros.*/
    int Bsize = DControllerMinPost.size();
    int count = 0;
    for (int i = 0; i < Bsize; i++)
    {
        if (i == 0) {
            CoarsePartition.push_back(DControllerMinPost[0]);
            count++;
            Bpartition_coarsePno[0] = count;           
        }
        else {
            for (int j = 0; j < count; j++) {
                if (DControllerMinPost[i] == CoarsePartition[j]) {
                    Bpartition_coarsePno[i] = j + 1;
                    break;
                }
            }
            if (Bpartition_coarsePno[i] == 0) {
                CoarsePartition.push_back(DControllerMinPost[i]);
                count++;
                Bpartition_coarsePno[i] = count;
            }
        }
    }
}

template<class state_type, class input_type, class F2, class F3, class F4>
void synthesis_invariance_mst(const std::vector<scots::abs_type>& domain_intersection,
    std::vector<bool>& still_safe,
    std::unique_ptr<bool[]>& input_list_status,
    scots::abs_type& Bsize,
    scots::UniformGrid& is,
    scots::UniformGrid& ss,
    scots::Abstraction<state_type, input_type>& abs,
    F2& system_post,
    F3& radius_post,
    F4& reach_set
    ) {
    /*domain_intersection = the safe set as state space grid indices*/
    scots::abs_type Dsize = domain_intersection.size();
    Bsize = Dsize;
  /* define function to check if the cell is in the safe set  */
    std::unordered_map<scots::abs_type, scots::abs_type> domain_intersection_xind_to_ind; // maps state-index to the index in the domain_intersection vector.
    for (scots::abs_type i = 0; i < Dsize; i++)
        domain_intersection_xind_to_ind[domain_intersection[i]] = i;
    auto safeset = [&still_safe, &domain_intersection_xind_to_ind](const scots::abs_type& idx) noexcept {
         if (domain_intersection_xind_to_ind.count(idx) > 0)
        {
            if (still_safe[domain_intersection_xind_to_ind[idx]])
                return true;
            else
                return false;
        }
        else
            return false;
    };
    
    /* compute winning domain (contains also valid inputs) */
    std::cout << "\nSynthesis for intersection: \n";

    scots::abs_type M = is.size();
        std::vector<std::unique_ptr<std::vector<scots::abs_type>>> post_list(M * Dsize);   // vector of unique pointers to vector of abs_type
        std::vector<int> n_v_inputs(Dsize, M);    // Number of valid inputs (initialized to M).
        for (scots::abs_type i = 0; i < Dsize * M; i++)
            input_list_status[i] = true;
        bool first_iteration = true;
        bool state_thrown_out = true; // Was any state thrown out of the safe set in the previous iteration.

        class s_grid_data_class {
        public:
            int dim;
            state_type eta, lower_left, upper_right;
            std::vector<scots::abs_type> NN;
        } s_grid_data_ob;
        s_grid_data_ob.dim = ss.get_dim();
        s_grid_data_ob.NN = ss.get_nn();
        for (int i = 0; i < s_grid_data_ob.dim; i++)
        {
            s_grid_data_ob.eta[i] = ss.get_eta()[i];
            s_grid_data_ob.lower_left[i] = ss.get_lower_left()[i];
            s_grid_data_ob.upper_right[i] = ss.get_upper_right()[i];
        }

        while (state_thrown_out)
        {
            state_thrown_out = false;
            for (scots::abs_type i = 0; i < Dsize; i++)
            {
                if (still_safe[i])
                {
                    for (scots::abs_type j = 0; j < M; j++)
                    {
                        if (first_iteration)
                        {
#if ReachSet == 0
                            /*post list may come out to be empty when the reachable set goes outside the state space grid.*/
                            abs.get_post_id(system_post, radius_post, s_grid_data_ob, domain_intersection[i], j, post_list[M * i + j]);
#else
                            //post_list[M * i + j] = &abs.get_post_id(reach_set, s_grid_data_ob, domain_intersection[i], j, temp_ob, "MST");
                            abs.get_post_id(reach_set, s_grid_data_ob, domain_intersection[i], j, post_list[M * i + j]);

#endif  // end ReachSet
                        }
                        if (input_list_status[M * i + j])
                        {
                            if (((post_list[M * i + j]) == nullptr)) // i.e. post went outside the state space grid.
                            {
                                input_list_status[M * i + j] = false;
                                n_v_inputs[i]--;
                                if (n_v_inputs[i] == 0)
                                {
                                    still_safe[i] = false;
                                    state_thrown_out = true;
                                    Bsize--;
                                }
                            }
                            else {
                                for (scots::abs_type k = 0; k < (post_list[M * i + j])->size(); k++)
                                {
                                    if (!safeset((*post_list[M * i + j])[k]))
                                    {
                                        /*The kth post went outside the safeset; invalidate the jth input. */
                                        post_list[M * i + j].reset();
                                        input_list_status[M * i + j] = false;
                                        n_v_inputs[i]--;
                                        if (n_v_inputs[i] == 0)
                                        {
                                            still_safe[i] = false;
                                            state_thrown_out = true;
                                            Bsize--;
                                        }
                                        break;
                                    }
                                }
                            }
                        }

                    }
                }
            }
            first_iteration = false;
        }
        //tt.toc();
        /*std::cout << "Winning domain size: " << Bsize << "\n";

        wfile_winDomain(still_safe, domain_intersection, input_list_status, ss, is, Example + "with_intersection_winDomain.scs", 1);
        std::cout << "\nWrite controller to " << Example << ".scs \n";
        wfile_controller_scs(still_safe, domain_intersection, input_list_status, ss, is, Example + ".scs");
        std::cout << "Done. \n";
        delete[] input_list_status;*/

        return;
}

template<class F1, class F2, class F3>
void wfile_winDomain(std::vector<bool>& still_safe, std::vector<scots::abs_type>& domain_intersection, F3& input_list_status, F1& ss, F2& is, std::string filename, int flag) {
    //flag: 1(write), 2(append)
    auto M = is.size();
    std::ofstream myfile;
    if(flag == 1)
        myfile.open(MATLABDIR + filename);
    else
        myfile.open(MATLABDIR + filename, std::ios::app);
    myfile << "#TYPE:WINNINGDOMAIN\n#SCOTS: i (state) j_0 ... j_n (valid inputs)\n#MATRIX:DATA\n#BEGIN:";
    myfile << ss.size() << " " << is.size() << "\n";
    for (scots::abs_type i = 0; i < still_safe.size(); i++)
        if (still_safe[i]) {
            myfile << domain_intersection[i] << " ";
            for (scots::abs_type j = 0; j < M; j++)
            {
                if (input_list_status[M * i + j])
                    myfile << j << " ";
            }
            myfile << "\n";
        }
    myfile << "#END\n";
    myfile.close();
}

template<class F1, class F2>
void wfile_winDomain_b(F1& ss, F2& is, const scots::abs_type& Bsize, const std::vector<scots::abs_type>& Bpartition, std::vector<std::vector<scots::abs_type>>& controllerMinPost, std::string filename, int flag) {
    //flag: 1(write), 2(append)
    std::ofstream myfile;
    if (flag == 1)
        myfile.open(MATLABDIR + filename);
    else
        myfile.open(MATLABDIR + filename, std::ios::app);
    if (!(myfile.is_open() )) {
        std::cout << "\n%%%%%%%%%%\nError: Failed to open file for write.\n%%%%%%%%%%";
        return ;
    }
    myfile << "#TYPE:WINNINGDOMAIN\n#SCOTS: i (state) j_0 ... j_n (valid inputs)\n#MATRIX:DATA\n#BEGIN:";
    myfile << ss.size() << " " << is.size() << "\n";
    for (scots::abs_type i = 0; i < Bsize; i++) {
        myfile << Bpartition[i] << " ";
        for (scots::abs_type j = 0; j < controllerMinPost[i].size(); j++)
            myfile << controllerMinPost[i][j] << " ";
        myfile << "\n";
    }
    myfile << "#END\n";
    myfile.close();
}

template<class F1, class F2>
void write_scots_header(F1& ss, F2& is, std::string filename) {
    int state_dim = ss.get_dim();
    int input_dim = is.get_dim();
    std::ofstream myfile;
    myfile.open(MATLABDIR+filename);
    myfile<< "#SCOTS:v0.2\n#TYPE:STATICCONTROLLER\n#SCOTS:STATE_SPACE\n#TYPE:UNIFORMGRID\n#MEMBER:DIM\n";
    myfile << state_dim <<"\n";
    myfile << "#VECTOR:ETA\n#BEGIN:" << state_dim << "\n";
    std::vector<double> ss_eta = ss.get_eta();
    for (auto eta_ith : ss_eta)
        myfile << eta_ith << "\n"; 
    myfile << "#END\n#VECTOR:LOWER_LEFT\n#BEGIN:" <<state_dim << "\n";
    auto ss_lower_left = ss.get_lower_left();
    for (auto aa : ss_lower_left)
        myfile << aa << "\n";
    myfile << "#END\n#VECTOR:UPPER_RIGHT\n#BEGIN:" << state_dim << "\n";
    auto ss_upper_right = ss.get_upper_right();
    for (auto aa : ss_upper_right)
        myfile << aa << "\n";
    myfile << "#END\n#SCOTS:INPUT_SPACE\n#TYPE:UNIFORMGRID\n#MEMBER:DIM\n" << input_dim ;
    myfile << "\n#VECTOR:ETA\n#BEGIN:" << input_dim<<"\n";
    auto is_eta = is.get_eta();
    for (auto eta_ith : is_eta)
        myfile << eta_ith << "\n";
    myfile << "#END\n#VECTOR:LOWER_LEFT\n#BEGIN:" << input_dim << "\n";
    auto is_lower_left = is.get_lower_left();
    for (auto aa : is_lower_left)
        myfile << aa << "\n";
    myfile << "#END\n#VECTOR:UPPER_RIGHT\n#BEGIN:" << input_dim << "\n";
    auto is_upper_right = is.get_upper_right();
    for (auto aa : is_upper_right)
        myfile << aa << "\n";
    myfile << "#END\n";
    //myfile << "#TYPE:WINNINGDOMAIN\n#SCOTS:i (state) j_0 ... j_n (valid inputs)\n#MATRIX:DATA\n#BEGIN:" << ss.size() << " " << is.size() << "\n";
    myfile.close();
}

template<class F1, class F2, class F3>
void wfile_controller_scs(std::vector<bool>& still_safe, std::vector<scots::abs_type>& domain_intersection, F3& input_list_status, F1& ss, F2& is, std::string filename) {
    write_scots_header(ss, is, filename);
    wfile_winDomain(still_safe, domain_intersection, input_list_status, ss, is, filename, 2);
}

template<class F1, class F2>
void wfile_controller_scs_b( F1& ss, F2& is, const scots::abs_type& Bsize, const std::vector<scots::abs_type>& Bpartition, std::vector<std::vector<scots::abs_type>>& controllerMinPost, std::string filename) {
    write_scots_header(ss, is, filename);
    wfile_winDomain_b(ss, is, Bsize, Bpartition, controllerMinPost, filename, 2);
}

void wfile(std::vector<std::vector<scots::abs_type>>& Xj, std::vector<scots::abs_type>& Bpartition) {
    unsigned long Bsize = Bpartition.size();
    std::unordered_map<scots::abs_type, int> Bpartition2indices;
    for (unsigned long i = 0; i < Bsize; i++)
        Bpartition2indices[Bpartition[i]] = i;  // indices start from 0.

    std::ofstream BiFile, BjFile;
    BiFile.open(MATLABDIR + "DataMatlab/S_dtCon_Bi.txt");
    BjFile.open(MATLABDIR + "DataMatlab/S_dtCon_Bj.txt");
    if (!(BiFile.is_open() && BjFile.is_open())) {
        std::cout << "\n%%%%%%%%%%\nError Writing Bi/Bj : Failed to open file for write.\n%%%%%%%%%%";
        return;
    }
    for (unsigned long i = 0; i < Bsize; i++) {
        for (unsigned long j = 0; j < Xj[i].size(); j++) {
            BiFile << i + 1 << " ";
            BjFile << Bpartition2indices[Xj[i][j]] + 1 << " ";           
        }
    }
    BiFile.close();
    BjFile.close();
    return;

}

template<class F>
void wfile(const std::vector<F>& V, const std::string str1, const int flag, const int plusOne, const std::string& Example){
    /*Writes the vector V to the file MATLABDIR+"DataMatlab/S_dtCon_"+str1+".txt".*/
        //flag: 1(write), 2(append)
        //plusOne: 1(plus 1), 0()
        std::ofstream myfile;
    
        if(flag==1)
            myfile.open(MATLABDIR+"DataMatlab/S_dtCon_"+str1+".txt");
        else
            myfile.open(MATLABDIR+"DataMatlab/S_dtCon_"+str1+".txt",std::ios::app);
        if (!myfile.is_open()) {
            std::cout << "\n%%%%%%%%%%\nError Writing " << str1 << " : Failed to open file for write.\n%%%%%%%%%%";
            return;
        }
        unsigned long Vsize = V.size();
        for(unsigned long i=0;i<Vsize;i++){
            if(plusOne)
                myfile << V[i]+1 << " ";
            else
                myfile << V[i] << " ";
        }
//        myfile << "];" << std::endl;
         myfile.close();
        return;
    }

template<class F1, class F2>
void wfile_DController(const std::vector<F1>& X, const std::vector<F2>& U, const std::string filename, const int flag, const int plusOne) {
    /*Writes the determinized controller to MATLABDIR/filenameDcontroller.txt*/
        //flag: 1(write), 2(append)
        //plusOne: 1(plus 1), 0()
    std::ofstream myfile;

    if (flag == 1)
        myfile.open(MATLABDIR + filename + "Dcontroller.txt");
    else
        myfile.open(MATLABDIR + filename + "Dcontroller.txt", std::ios::app);
    if (!myfile.is_open()) {
        std::cout << "\n%%%%%%%%%%\nError Writing Determinized controller : Failed to open file for write.\n%%%%%%%%%%";
        return;
    }
    unsigned long dlength = X.size();
    for (unsigned long i = 0; i < dlength; i++) {
        if (plusOne)
            myfile << X[i] + 1 << " " << U[i] +1 << "\n";
        else
            myfile << X[i] << " " << U[i] << "\n";
    }
     myfile.close();
    return;
}

template<class F>
void wfile(const std::vector<std::vector<F>>& V, const int flag, const int plusOne, const std::string& Example){
    //flag: 1(write), 2(append)
    //plusOne: 1(plus 1), 0()
    std::ofstream BiFile, BjFile;
    if(flag==1){
        BiFile.open(MATLABDIR+"DataMatlab/S_dtCon_Bi.txt");
        BjFile.open(MATLABDIR+"DataMatlab/S_dtCon_Bj.txt");
    }
    else{
        BiFile.open(MATLABDIR+"DataMatlab/S_dtCon_Bi.txt",std::ios::app);
        BjFile.open(MATLABDIR+"DataMatlab/S_dtCon_Bj.txt",std::ios::app);
    }
    if (!(BiFile.is_open() && BjFile.is_open())) {
        std::cout << "\n%%%%%%%%%%\nError Writing Bi/Bj : Failed to open file for write.\n%%%%%%%%%%";
        return;
    }
    unsigned long Vsize = V.size();
    for(unsigned long i=0; i<Vsize; i++){
        for(unsigned long j=0; j<V[i].size(); j++){
            if(plusOne){
                BiFile << i+1 << " ";
                BjFile << V[i][j] +1 << " ";
            }
            else{
                BiFile << i << " ";
                BjFile << V[i][j]  << " ";
            }
          }
    }
    BiFile.close();
    BjFile.close();
    return;

}

template<class F>
void wfile(const std::vector<std::vector<F>>& Xj, const std::vector<F>& Bpartition, const int flag, const int plusOne, const std::string& Example) {
    //flag: 1(write), 2(append)
    //plusOne: 1(plus 1), 0()
    std::ofstream XiFile, XjFile;
    if (flag == 1) {
        XiFile.open(MATLABDIR + "DataMatlab/Xi.txt");
        XjFile.open(MATLABDIR + "DataMatlab/Xj.txt");
    }
    else {
        XiFile.open(MATLABDIR + "DataMatlab/Xi.txt", std::ios::app);
        XjFile.open(MATLABDIR + "DataMatlab/Xj.txt", std::ios::app);
    }
    if (!(XiFile.is_open() && XjFile.is_open())) {
        std::cout << "\n%%%%%%%%%%\nError Writing Xi/Xj : Failed to open file for write.\n%%%%%%%%%%";
        return;
    }
    unsigned long Xisize = Bpartition.size();
    for (unsigned long i = 0; i < Xisize; i++) {
        for (unsigned long j = 0; j < Xj[i].size(); j++) {
            if (plusOne) {
                XiFile << Bpartition[i] + 1 << " ";
                XjFile << Xj[i][j] + 1 << " ";
            }
            else {
                XiFile << Bpartition[i] << " ";
                XjFile << Xj[i][j] << " ";
            }
        }
    }
    XiFile.close();
    XjFile.close();
    return;

}


void wfile(const double &T, const std::string str1, const int flag, const int plusOne, const std::string& Example){
        //flag: 1(write), 2(append)
        //plusOne: 1(plus 1), 0()
        std::ofstream myfile;
    
        if(flag==1)
            // myfile.open("/Users/mst/GoogleDrive/Currently_working_on_these/Matlab/S_dtCon_"+str1+".txt");
            myfile.open(MATLABDIR+"DataMatlab/S_dtCon_"+str1+".txt");
        else
            // myfile.open("/Users/mst/GoogleDrive/Currently_working_on_these/Matlab/S_dtCon_"+str1+".txt",std::ios::app);
            myfile.open(MATLABDIR+"DataMatlab/S_dtCon_"+str1+".txt",std::ios::app);
        if (!myfile.is_open()) {
            std::cout << "\n%%%%%%%%%%\nError Writing " << str1 << " : Failed to open file for write.\n%%%%%%%%%%";
            return;
        }

    if(plusOne)
        myfile << T+1 << " ";
    else
        myfile << T << " ";
        
//        myfile << "];" << std::endl;
    myfile.close();
    return;
    }

template<class F, class F2>
void wfile_pFacesSafeset(F& ss, F2& eta, std::vector<scots::abs_type>& D) {
    /*saves the domain as lb,ub of cells for pFaces safeset.*/
    int num = D.size();
    F2 cell_centre;
    const int state_dim = cell_centre.size();
    std::ofstream myfile;
    myfile.open(MATLABDIR + "pFaces_Safeset.txt");
    //myfile.precision(16);
    myfile << "safe{\n    count=\"" << num << "\";\n";
    for (int i = 0; i < num; i++)
    {
        myfile << "    s" << i + 1 << "{type=\"rectangle\"; h=\"";
        ss.itox(D[i], cell_centre);
        for (int j = 0; j < state_dim; j++)
        {
            myfile << "{";
            myfile << cell_centre[j] - 0.5 * eta[j] <<"," << cell_centre[j] + 0.5 * eta[j] <<"}";
            if (j != state_dim - 1)
                myfile << ",";
        }
        myfile << "\";}\n";
    }
    myfile << "}";
    myfile.close();
    return;
}

template<class F>
void wfile(const std::vector<F>& B, std::string filename){
    // std::ofstream myfile;
    // myfile.open(filename);
    // for(unsigned long i=0;i<B.size(); i++)
        // myfile << B[i] << "\n";
    return;
}

template<class F>
bool read_winDomain_mst(std::string filename, std::vector<std::vector<F>>& ivec, std::vector<F>& domain) {
    F val;   
    std::string line;
    std::ifstream infile(filename);
    if (infile.is_open())
    {
        for (int i = 0; i < 4; i++)
            std::getline(infile, line);
        while (std::getline(infile, line))
        {
            if (line == "#END")
                return true;    // end of the scs file.
            std::vector<F> ilist;   // input list
            std::istringstream iss(line);
            iss >> val;
            domain.push_back(val);
            while (iss >> val)
                ilist.push_back(val);
            ivec.push_back(ilist);
        }
        return true;
    }
    else
        return false;
}

template<class F>
bool rfile(std::string filename, std::vector<F> &vec){
	std::string line;
	std::ifstream infile(filename);
	std::getline(infile, line);
	std::istringstream iss(line);
	F val;
	while(iss >> val)
		vec.push_back(val);
	return true;
}



#endif /* Helper_MST_h */
