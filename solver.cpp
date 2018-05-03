/* maxflow.c */

/* Written by Andrew Makhorin <mao@gnu.org>, October 2015. */

#include <math.h>


#include <glpk.h>

#include <queue>
#include <set>

#include <time.h>
#include <cstdlib>

#include <vector>
#include <list>
#include <utility>
#include <cstdio>



#include <iostream>
#include <algorithm>

#include <cmath>
#include <sstream>
#include <fstream>
#include <string>

#include <graphics.h>

using namespace std;

const int clauses_num = 3;
const int variables_num = 5;
int active_variables [clauses_num][3] ;

int * ia = new int[50];
int * ja = new int[50];
double * ar = new double[50];



void initialize_with_input(const char * file, int * ia, int *ja, double * ar){
  

  int j;
  int step_count =1;

  int rows = 0;
  ifstream input(file);
  if (!input.is_open())
    cout << "err opening file " <<endl;
  else{
    string line;
    while (getline(input, line))
    {
        stringstream curr_line(line);
        int n; 
        int var = 0;
        while(curr_line >> n) 
          {
              active_variables[rows][var] = n;
              var++;

          }

              for (j =1; j<=variables_num; j++){


                ia[step_count] = rows+1;

                ja[step_count] = j;

                //set active variables to have coefficients 1
                if ( j==active_variables[rows][0] || j==active_variables[rows][1] || j==active_variables[rows][2]){
                      ar[step_count] = 1;
                }
                else
                      ar[step_count] = 0;

                step_count+=1;

              }
        // process pair (a,b)
          rows++;
    }

  }
}

void initialize_with_fixed_input(int * ia, int *ja, double * ar, int flag){
  

  int j;
  int step_count =1;

  int rows = 0;

  ifstream input;

cout<< "FLAG:  " << flag <<endl;
  if (flag == 1)
    input.open("input.txt");
  else if (flag == 2)
    input.open("input2.txt");
  else {
    cout << "error file" <<endl;
    return;
  }
  if (!input.is_open())
    cout << "err opening file " <<endl;
  else{
    string line;
    while (getline(input, line))
    {
        stringstream curr_line(line);
        int n; 
        int var = 0;
        while(curr_line >> n) 
          {
              active_variables[rows][var] = n;
              var++;

          }

              for (j =1; j<=variables_num; j++){


                ia[step_count] = rows+1;

                ja[step_count] = j;

                //set active variables to have coefficients 1
                if ( j==active_variables[rows][0] || j==active_variables[rows][1] || j==active_variables[rows][2]){
                      ar[step_count] = 1;
                }
                else
                      ar[step_count] = 0;

                step_count+=1;

              }
        // process pair (a,b)
          rows++;
    }

    input.close();
  }
}



void initialize_random(int *& ia, int *&ja, double *& ar){

  vector<int> pool1;
  vector<int> pool2;
  vector<int> pool3;


  //set up constraints (0 <= clause <= 1)
  //set up how we create clauses
  for (int i =1; i<=variables_num; i++){

        pool1.push_back(i);
        pool2.push_back(i);
        pool3.push_back(i);
  }


  int step_count =1;



  //randomly generate clauses for our problem instance
  random_shuffle(pool1.begin(), pool1.end());
  random_shuffle(pool2.begin(), pool2.end());
  random_shuffle(pool3.begin(), pool3.end());


  for (int i =1; i<=clauses_num; i++){


        int r1 = pool1.back();

        pool1.pop_back();

        int r2 = pool2.back();

        //avoid same variable being picked twice for same clause
        if (r2 == r1){ 
              r2 = pool2.at(pool2.size()-2);
              pool2.erase(pool2.end()-2);
        }
        else{
              pool2.pop_back();
        }

        int r3 = pool3.back();


        if (r3 == r1 || r3 == r2){ 
                    
              r3 = pool3.at(pool3.size()-2);


              if (r3== r1 || r3 == r2){


                    r3 = pool3.at(pool3.size()-3);
                    pool3.erase(pool3.end()-3);

              }

              else{
                    pool3.erase(pool3.end()-2);

              }

        }
        else{
              pool3.pop_back();
        }


        //keep track of what variables are active in a single clause
        active_variables[i-1][0] = r1;
        active_variables[i-1][1] = r2;
        active_variables[i-1][2] = r3;



        for (int j =1; j<=variables_num; j++){


              ia[step_count] = i;

              ja[step_count] = j;

              //set active variables to have coefficients 1
              if ( j==r1 || j== r2 || j==r3){
                    ar[step_count] = 1;
              }
              else
                    ar[step_count] = 0;

              step_count+=1;

        }


  }
}

void setup_empty_lp(glp_prob *& lp){

//create problem through glpk library, each row is clause,
//each column is one variable
//we minimize answer
lp = glp_create_prob();

glp_set_prob_name(lp, "sample");

glp_set_obj_dir(lp, GLP_MIN);

glp_add_rows(lp, clauses_num);



for (int i =1; i<=clauses_num; i++){
      char temp_clause_name [20]= "clause";
      sprintf(temp_clause_name+6, "%d", i);

      glp_set_row_name(lp, i, temp_clause_name);

      glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0);
}



      glp_add_cols(lp, variables_num);

      //in minimization problem, we choose to minimize
      //simply the expression x1 
      // doesn't really matter, it will all be correct in end, 
      // just need to determine correct value of one value (or if imppssible)
      glp_set_obj_coef(lp, 1, 1.0);



  for (int i =1; i<=variables_num; i++){

    char temp_var_name [12]= "x";
    sprintf(temp_var_name+1, "%d", i);

    glp_set_col_name(lp, i, temp_var_name);

    glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);

    glp_set_obj_coef(lp, i, 0.0);

  }

}

void create_graph(double * solutions, vector< list< pair<int, double> > > &  adjacencyList){
        for (int i = 1; i <= variables_num + clauses_num ; i++) {
         
            //add pairs of active variables, current solution of variables to ith clause
            if (i <= clauses_num){
                  for (int j = 0; j<3; j++){

                        double temp = solutions[ active_variables[i-1][j]];

                        adjacencyList[i].push_back(make_pair( active_variables[i-1][j]+clauses_num, temp ));
                  }
            }
            
            //add pairs of clauses ith variable is in, current solution of the variable 
            if ( i > clauses_num){

                  for (int k = 0; k < clauses_num; k++){

                        if ( active_variables[k][0] == i - clauses_num || active_variables[k][1] == i - clauses_num || active_variables[k][2] == i - clauses_num  ){
                              
                              double temp =  solutions[i - clauses_num];

                              adjacencyList[i].push_back(make_pair( k+1, temp ));
                        }

                        
                  }  
            }
      }
}


void draw_graph(double * solutions, long long int **  new_solids, int rows_num, int clauses_num){
          int gd = DETECT, gm;
        initgraph(&gd, &gm, NULL);

        //initwindow(400, 300);
 
          for (int i=1; i<= clauses_num; i++){

            for (int j=1; j<=variables_num; j++){


              if (active_variables[i-1][0] == j || active_variables[i-1][1] == j || active_variables[i-1][2] == j){
                  int nearest = (int) (solutions[j]* 25);  
                  setcolor(i);
                  circle(50*j, 50*i, nearest);

              }



            }

          }
                  setcolor(15);

                  line(0, 175, 400, 175);


          for (int i=1; i<= rows_num-clauses_num; i++){



            for (int j=1; j<=variables_num+1; j++){

              int nearest = (int) (new_solids[i][j]* 5);  
              setcolor(10);
              circle(50*j, 50*(i+clauses_num), nearest);

            }

          }

        //delay(10000);
        //Wait for a key press
        int in = 0;

        while (in == 0) {
            in = getchar();
        }

        closegraph();


}



void work_loop(double * solutions, glp_prob *& lp, int rows_num, int clauses_num, long long int ** new_solids, int flag){

      //use this later in work look
      std::queue<int> graph_queue;
      std::queue<int> variables_queue;
      std::queue<int> clauses_queue;

      set<short int> border_set;

      set<short int> border2_set;

  double z;
bool flag_border = true;

//main work loop
while (flag_border){

      border_set.clear();
      border2_set.clear();

      int i,j, k;

      //get current solutions and print them out to console
      for (i =1; i<=variables_num; i++){
           solutions[i] = glp_get_col_prim(lp, i);
           //cout << "soln " << solutions[i] << endl;
      }

      for (i=1; i<= variables_num; i++){

           printf("| %d = %g | ",
            i, solutions[i]);

      }


      //represent our current problem as graph through adjacency list
      vector< list< pair<int, double> > > adjacencyList(variables_num + clauses_num + 1);
      create_graph(solutions, adjacencyList);


      //set up our temp solutions array
      int fractional_solutions [variables_num+1] ;

      for (i=0; i< variables_num +1; i++)
      {      
        fractional_solutions[i]= 0;
      }


      for (i = 1; i <= clauses_num; i++) {
              //printf("adjacencyList[%d] ", i);
               
            list< pair<int, double> >::iterator itr = adjacencyList[i].begin();

            while (itr != adjacencyList[i].end()) {
                  //because of rounding issues, make nonzero into 1's
                  if ( abs((*itr).second - 1.0) > 0.0000001 && abs((*itr).second) > 0.0000001 ){

                        fractional_solutions[(*itr).first-clauses_num] = 1;

                  }
                  
                  ++itr;
            }
        

      }


      //initialzed all vertices as unvisited
      int visited [variables_num+ clauses_num + 1];

      for (i=0; i< variables_num+clauses_num+1; i++){
            visited[i]= 0;
      }

      int border [variables_num+ 1];

      //border elements are the variables that have value 0 on current 
      // iteration but in at least one clause, the other two variables happen
      // to be nonzero 
      // (aka a 1,0,0 clause is good  )
      // 0.5,0.5,0 clause is bad, we add the last variable to the border set 
      for (i=0; i< variables_num+1; i++){
            border[i]= 0;
      }




      int border_counter = 0;

    //std::queue<int>().swap(clauses_queue);
    //std::queue<int>().swap(variables_queue);
    //std::queue<int>().swap(graph_queue);



      bool flag1 =false;
      bool flag2 = false;

      short int temp_border;


      //initialzed all vertices as unvisited
      for (i = 1; i <= variables_num; i++) {

            for (j=0; j< variables_num+clauses_num+1; j++){
                   visited[j]= 0;
            }
            //if a solution is zero, it's "good", or solved, so we don't need to add constraints to it 
            //else we do below




            if (fractional_solutions[i] ==1 ){

                  visited[i] = 1;

                  variables_queue.push(clauses_num+i);

                  //bfs over all vertices to add clauses with non-zero/ non-one variables
                  while (variables_queue.size()>0){



                        while (variables_queue.size()>0){


                              //cout << "vq size:  " << variables_queue.size() <<endl;
                              //if (!variables_queue.empty())
                              //  cout << "vq val:  " << variables_queue.front() <<endl;
                              int curr_variable = variables_queue.front();

                              if (!variables_queue.empty()){
                                variables_queue.pop();
                              }
                              
                        
                              list< pair<int, double> >::iterator itr = adjacencyList[curr_variable].begin();

                              //find unvisited nonzero value variable,
                              // we want to check all of variables clauses to 
                              // see if any border elements exist
                              while (itr != adjacencyList[curr_variable].end()) {
                                    if ( visited[(*itr).first] != 1 &&  abs((*itr).second-1.0) > 0.00001 && abs((*itr).second) > 0.00001 ){

                                          int temp = (*itr).first;
                                          clauses_queue.push(temp);
                                          visited[(*itr).first] = 1 ;

                                    }
                                    ++itr;
                              }

                              /*
                              cout << "vq size2:  " << variables_queue.size() <<endl;
                              if (!variables_queue.empty())
                                cout << "vq val2:  " << variables_queue.front() <<endl;
                              */
                        }


                        //go through all clauses with nonzero/ nonone solutions
                        while (clauses_queue.size()>0){


                              int curr_clause;

                              /*
                              cout << "cq size:  " << clauses_queue.size() <<endl;
                              if (!clauses_queue.empty())
                                cout << "cq val:  " << clauses_queue.front() <<endl;
                              */
                              curr_clause = clauses_queue.front();
  
                              if (!clauses_queue.empty()){
                                clauses_queue.pop();
                              }


                              list< pair<int, double> >::iterator itr = adjacencyList[curr_clause].begin();

                              flag1 = false;
                              flag2 = false;
                              //int nonzero_triple[3];

                              //go through all variables in "interesting"  clause
                              // aka potentially border elements
                              while (itr != adjacencyList[curr_clause].end()) {

                                    //there is zero element in clause, set that variable as border for lp
                                    if ( abs((*itr).second) < 0.00001 ){

                                          flag1 = true;
                                          temp_border = (short) (*itr).first - clauses_num;
                                    }


                                    if (  abs((*itr).second-1.0) > 0.00001 && abs((*itr).second) > 0.00001  ){

                                          flag2 = true;
                                    }


                                    if (visited[(*itr).first] != 1 &&  abs((*itr).second-1.0) > 0.00001 && abs((*itr).second) > 0.00001  ){
                                          variables_queue.push((*itr).first);
                                          visited[(*itr).first] = 1;
                                          //cout << "INSIDE" <<endl;

                                    }

                                    /*
                                  cout << "vq size3:  " << variables_queue.size() <<endl;
                                  if (!variables_queue.empty())
                                    cout << "vq val3:  " << variables_queue.front() <<endl;
                                    
                                  cout << "BS size:  " << border_set.size() <<endl;
                                  if (!border_set.empty())
                                    cout << "BS val:  " << *(border_set.begin()) <<endl;
                                    */

                                    //we have found a border clause
                                    // add it to the border set data structure
                                    if (flag1 && flag2){
                                        
                                        //cout << "INSERT BS" <<endl;
                                        //cout << "vq sizeBS1:  " << variables_queue.size() <<endl;
                                        //cout << "cq sizeBS1:  " << clauses_queue.size() <<endl;
                                        //cout << "BS sizeBS1:  " << border_set.size() <<endl;

                                        border[border_counter] = temp_border; 
                                        border_set.insert(temp_border);

                                        border_counter+=1;

                                    }
                                    /*
                                  cout << "BS size2:  " << border_set.size() <<endl;
                                  if (!border_set.empty())
                                    cout << "BS val2:  " << *(border_set.begin()) <<endl;

                                  cout << "vq size4:  " << variables_queue.size() <<endl;
                                  if (!variables_queue.empty())
                                    cout << "vq val4:  " << variables_queue.front() <<endl;
                                    */
                  
                                    ++itr;
                              }


                        }
                  }


            }

      }


      /*
      set <short int> :: iterator itr;

      for (itr = border_set.begin(); itr != border_set.end(); ++itr) {
            //cout<< " border " << (*itr);
      }
  */
      rows_num++;

      //no border means we have found a all 1s and all 0s
      //solution, aka we found valid solution to problem and 
      //we are done
      if (border_set.size() == 0){


            cout << "NO BORDER" <<endl;
            flag_border = false;
            break;

      }

      //otherwise we add the border elements as another constraint

      //add constraint
      glp_add_rows(lp, 1);


      char temp_clause_name [12]= "solid";
      sprintf(temp_clause_name+5, "%d", rows_num- clauses_num );

      glp_set_row_name(lp, rows_num, temp_clause_name);



      long long int current_border_number = 0;

      int border_cumulative = 0;

      long double current_border_sum = 0.0;

      unsigned long long int current_border_sum2 = 0;



      int current_solid [variables_num];

      int current_border_length = 0;




      set<int> exact_border_set;

      long long int border_element_frequency[variables_num+1];

      //this array will check how much border variable came up before
      //use that to increment new border (so we make progress)
      for (i=1; i<= variables_num; i++){
            border_element_frequency[i] = 0;
      }




      for (i=1; i<= rows_num-clauses_num-1; i++){

            current_border_length = 0;
            current_border_sum = 0.0;

            current_border_sum = 0;


            //check if current solution exactly satisfies one of previous 
            // border constraints

            for (j=1; j<=variables_num; j++){

                  if (  new_solids[i][j] != 0 ){

                        
                        current_solid[j] = 1;

                        //]]cout << "pro  " <<  solutions[j] << " : " << new_solids[i][j] <<endl;

                        long double nearest = roundf(solutions[j]*new_solids[i][j] * 10000) / 10000;  

                        current_border_sum +=  nearest;


                   }

            }


            //if it does, we add that equation to the new constraint we are making
            //this is to make faster progress (hopefully)
            if (abs(current_border_sum - new_solids[i][variables_num+1]) < 0.0001){

                  current_border_number+=1;

                  //we add the previous bounds, and the coefficients of the variables
                  border_cumulative += new_solids[i][variables_num+1];

                  for(k=1; k<= variables_num; k++){


                        if (  new_solids[i][k] != 0 ){

                              exact_border_set.insert(k);


                              border_element_frequency[k] += new_solids[i][k];

                              //cout << " INSERT FREQ AT " << k << endl;

                         }

                  }
            }
            
      }


      //increment current border elements, if not appear before
      //make them 1
      for (i=1; i<= variables_num; i++){

            if ( (border_set.find(i)!= border_set.end()) && (exact_border_set.find(i)!= exact_border_set.end()) ){
                 new_solids[rows_num-clauses_num][i] = border_element_frequency[i] + 1;

             }

             else if ( exact_border_set.find(i) != exact_border_set.end() ){
                 new_solids[rows_num-clauses_num][i] = border_element_frequency[i];

             }

             else if ( border_set.find(i)!= border_set.end()  ){
                 new_solids[rows_num-clauses_num][i] =1;

             }


            else {
                 new_solids[rows_num-clauses_num][i] =0;
            }


      }

      new_solids[rows_num-clauses_num][variables_num+1] = 1+border_cumulative;



      if (flag==0)
        draw_graph(solutions, new_solids, rows_num, clauses_num);

      glp_set_row_bnds(lp, rows_num, GLP_LO, 1+border_cumulative, 1+border_cumulative);



      //set up new lp problem
      for (i =rows_num; i<=rows_num; i++){

            for (j =1; j<=variables_num; j++){



                  ia[(rows_num-1)*variables_num+j] = i;

                  ja[(rows_num-1)*variables_num+j] = j;

                  if ( (border_set.find(j) != border_set.end() ) || ( exact_border_set.find(j) != exact_border_set.end() )) {


                    if ( (border_set.find(j) != border_set.end() ) && ( exact_border_set.find(j) != exact_border_set.end() )) {
                        //cout << "BOTH  " <<  border_element_frequency[j] << endl;
                         ar[(rows_num-1)*variables_num+j] = border_element_frequency[j]+1;
                    }

                  else if ( (border_set.find(j) != border_set.end() ) ) {
                         ar[(rows_num-1)*variables_num+j] = 1;
                    }

                    else {
                         ar[(rows_num-1)*variables_num+j] = border_element_frequency[j];

                    }


                       
                  }

                  else 
                        ar[(rows_num-1)*variables_num+j] = 0;

            }


      }


      //get solutons back
      glp_load_matrix(lp, variables_num*rows_num, ia, ja, ar);


      //new instance of problem (with additional constraint added )
      glp_write_lp(lp, NULL, "problem2.txt");


      glp_simplex(lp, NULL);

      if (glp_get_status(lp) == GLP_NOFEAS){
            break;
      }

      z = glp_get_obj_val(lp);


}

}


int * solve(int argc, char *argv[], int flag = 0) {


int i;
int j;
int k;

glp_prob *lp;

double z, x1, x2, x3;


//used to fill up array to pass in to lp library


//number of constraints we fan add
const int new_solids_length = 50;





long long int ** new_solids = new long long int*[new_solids_length];
for(i = 0; i < new_solids_length; i++)
    new_solids[i] = new long long int[new_solids_length];


for (int i =0; i< new_solids_length; i++){
      for (int j =0; j< new_solids_length; j++){
            new_solids[i][j] = 0;
      }
}





srand(time(NULL));   // should only be called once




int rows_num = clauses_num;


setup_empty_lp(lp);

  if (flag == 1){
    initialize_with_fixed_input( ia, ja, ar, 1);
  }

  else if (flag == 2){
    initialize_with_fixed_input( ia, ja, ar, 2);

  }


  else if (argc==1){
    initialize_random( ia, ja, ar);
  }


  else if (argc==2){
    initialize_with_input(argv[1], ia, ja, ar);
  }





FILE* fp = stdout;

//load the lp library with our problem instance
glp_load_matrix(lp, variables_num*clauses_num, ia, ja, ar);

//original instance of problem
glp_write_lp(lp, NULL, "problem.txt");

glp_simplex(lp, NULL);

z = glp_get_obj_val(lp);



//set solutions array
double * solutions  = new double[variables_num+1];
      for (i =0; i<=variables_num; i++){
           solutions[i] = 0.0;

      }


work_loop(solutions, lp, rows_num, clauses_num, new_solids, flag);


//prepare output and clean up memory
int * res = new int[variables_num];

for(i = 0; i < new_solids_length; i++)
    delete[] new_solids[i];
delete[] new_solids;

delete[] ia;
delete[] ja;

delete[] ar;

for (i =0; i<variables_num; i++)
  res[i] =solutions[i+1];
delete[] solutions;
glp_delete_prob(lp);
glp_free_env();

  return res;

}

