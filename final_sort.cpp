#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include <vector>
#include <sstream>

#include <limits>
#include <algorithm>
#include <utility>
using namespace std;
//mpiCC final_sort.cpp -o sort2d2
//mpirun inputfile outputfile rows columns
// mpirun sort2d2 Lab4_tests/input1 Lab4_tests/output1 100 100

bool compare_as_floats (float i,float j)
{
  return ((i)<(j));
}



int main(int argc,char** argv){

	string inputfile=argv[1];
	string outputfile=argv[2];
	int maxrows=atoi(argv[3]);
	int maxcols=atoi(argv[4]);


	vector<int> row_v_key1;
	vector<float> row_v_key2;
	vector<int> row_v_index;
	vector<int> row_v_countcum;

	vector<int> col_v_key1;
	vector<float> col_v_key2;
	vector<int> col_v_index;
	vector<int> col_v_countcum;
	

	MPI_Init(&argc,&argv);
	int world_size;
	int world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	
	if (world_rank==0){
		std::vector<std::vector<std::pair<int,int> > > v_key1(maxrows, std::vector<std::pair<int,int> > (0));
		std::vector<std::vector<std::pair<int,float> > > v_key2(maxrows, std::vector<std::pair<int,float> > (0));
		FILE *fp;
	   	fp = fopen(inputfile.c_str(), "rb");
		if(fp == NULL){
	        	printf("Error opening file\n");
	        	exit(1);
	    	}

		int row;
		int col;
		int key1;
		float key2;
		int count=0;
	    	while( fread(&row,sizeof(int), 1, fp) == 1 )
	    	{
			fread(&col, sizeof(int), 1, fp);
			fread(&key1,sizeof(int), 1, fp);
			fread(&key2,sizeof(float), 1, fp);
			v_key1[row].push_back(make_pair(col,key1));
			v_key2[row].push_back(make_pair(col,key2));
			count+=1;
	    	}
	    	fclose(fp);
		
		int cumulative=0;
		for (int i=0;i<maxrows;i++){
			vector<std::pair<int,int> > valkey1 = v_key1[i]; 
			vector<std::pair<int,float> > valkey2 = v_key2[i];

 			stable_sort(valkey1.begin(),valkey1.end());
			stable_sort(valkey2.begin(),valkey2.end());
			for (int j=0;j<valkey1.size();j++){
				row_v_key1.push_back(valkey1[j].second);
				row_v_key2.push_back(valkey2[j].second);
				row_v_index.push_back(valkey1[j].first);
			}

			row_v_countcum.push_back(cumulative);
			cumulative += valkey1.size();
		}
		row_v_countcum.push_back(cumulative);
		v_key1.clear();
		v_key2.clear();
		
	}

	for (int loop=0;loop<4;loop++){
		// SORT ROW WISE
		int size= maxrows/world_size;
		if (world_rank==0){
			for (int i=1;i<world_size;i++){
				int start = row_v_countcum[(i*size)];
				int end = row_v_countcum[(i+1)*size];
				int length = end-start;
	 			vector<int> all_vals_key1;
	 			vector<float> all_vals_key2;
				int counter=0;
				for (int j=start;j<end;j++){
					all_vals_key1.push_back(row_v_key1[j]);
					all_vals_key2.push_back(row_v_key2[j]);
				}

				MPI_Send(&length,1,MPI_INT,i,0,MPI_COMM_WORLD);
	 			MPI_Send(&row_v_countcum[i*size],size+1,MPI_INT,i,1,MPI_COMM_WORLD);
	 			MPI_Send(&all_vals_key1[0],length,MPI_INT,i,2,MPI_COMM_WORLD);	
	 			MPI_Send(&all_vals_key2[0],length,MPI_FLOAT,i,3,MPI_COMM_WORLD);

			}
		}

		int vectlength=0;
		if (world_rank!=0){
		  	MPI_Recv(&vectlength,1, MPI_INT, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		vector < int > v1(size+1,0);
		vector < int > v2(vectlength,0);
		vector < float > v3(vectlength,0);
		if (world_rank!=0){
			MPI_Recv(&v1[0],size+1, MPI_INT, 0,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&v2[0],vectlength, MPI_INT, 0,2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  	MPI_Recv(&v3[0],vectlength, MPI_FLOAT, 0,3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			vector< std::pair <float,int> > vect;
			for (int i=0;i<vectlength;i++){
				vect.push_back(make_pair(v3[i],v2[i]));
			}
			
			int indexstart = v1[0];
			for (int i=0;i<size;i++){
				vector< std::pair <float,int> > keys;
				int start = v1[i];
				int end = v1[i+1];
				for (int j=start-indexstart;j<end-indexstart;j++){
					keys.push_back(vect[j]);
				}
				stable_sort(keys.begin(),keys.end(),[](const std::pair<float, int>& p1, const std::pair<float, int>& p2) { return p1.first < p2.first;});	
				for (int j=0;j<keys.size();j++){
					vect[start-indexstart+j] = keys[j];
				}	
			}
			for (int i=0;i<vect.size();i++){
				v2[i]=vect[i].second;
				v3[i]=vect[i].first;
			}
			MPI_Send(&vectlength,1,MPI_INT,0,0,MPI_COMM_WORLD);
			MPI_Send(&v2[0],vectlength,MPI_INT,0,1,MPI_COMM_WORLD);	
	 		MPI_Send(&v3[0],vectlength,MPI_FLOAT,0,2,MPI_COMM_WORLD);
		}


		if (world_rank==0){

			for (int i=1;i<world_size;i++){
				int vectlength;
		  		MPI_Recv(&vectlength,1, MPI_INT, i,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				vector<int> all_vals_key1(vectlength,0);
	 			vector<float> all_vals_key2(vectlength,0);
				MPI_Recv(&all_vals_key1[0],vectlength, MPI_INT, i,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  		MPI_Recv(&all_vals_key2[0],vectlength, MPI_FLOAT, i,2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				int startpoint= row_v_countcum[i*size];				
				for (int j=0;j<vectlength;j++){
					row_v_key1[startpoint+j]=all_vals_key1[j];
					row_v_key2[startpoint+j]=all_vals_key2[j];
				}		
			}
			
			
			
			// SORT ITS OWN PROCESS DATA
			//INITIAL ROWS
			for (int i=0;i<size;i++){
				int start = row_v_countcum[i];
				int end = row_v_countcum[(i+1)];
				int length = end-start;
				vector< std::pair <float,int> > vect;
				for (int j=start;j<end;j++){
					vect.push_back(make_pair(row_v_key2[j],row_v_key1[j]));
				}
				stable_sort(vect.begin(),vect.end(),[](const std::pair<float, int>& p1, const std::pair<float, int>& p2) { return p1.first < p2.first;});	
				for (int j=0;j<vect.size();j++){
					row_v_key1[start+j]=vect[j].second;
					row_v_key2[start+j]=vect[j].first;
				}
			}

			//END ROWS
			for (int i=world_size*size;i<maxrows;i++){
				int start = row_v_countcum[i];
				int end = row_v_countcum[(i+1)];
				int length = end-start;
				vector< std::pair <float,int> > vect;
				for (int j=start;j<end;j++){
					vect.push_back(make_pair(row_v_key2[j],row_v_key1[j]));
				}
				stable_sort(vect.begin(),vect.end(),[](const std::pair<float, int>& p1, const std::pair<float, int>& p2) { return p1.first < p2.first;});	
				for (int j=0;j<vect.size();j++){
					row_v_key1[start+j]=vect[j].second;
					row_v_key2[start+j]=vect[j].first;
				}
			}

			int cumulative=0;

			for (int i=0;i<maxcols;i++){
				int temp_cumulative=cumulative;
				for (int j=0;j<maxrows;j++){
					int start = row_v_countcum[j];
					int end = row_v_countcum[(j+1)];
					int length = end-start;
					for (int k=0;k<length;k++){
						if (row_v_index[start+k]==i){
							col_v_key1.push_back(row_v_key1[start+k]);
							col_v_key2.push_back(row_v_key2[start+k]);
							col_v_index.push_back(j);
							cumulative+=1;
						}
					}
				
				} 
				col_v_countcum.push_back(temp_cumulative);
				
			}
			col_v_countcum.push_back(cumulative);
			row_v_index.clear();
			row_v_key1.clear();
			row_v_key2.clear();
			//row_v_countcum.clear();
		}

		

		// NOW DO COLUMN WISE SORTING

		size= maxcols/world_size;
		if (world_rank==0){
			for (int i=1;i<world_size;i++){
				int start = col_v_countcum[(i*size)];
				int end = col_v_countcum[(i+1)*size];
				int length = end-start;
	 			vector<int> all_vals_key1;
	 			vector<float> all_vals_key2;
				for (int j=start;j<end;j++){
					all_vals_key1.push_back(col_v_key1[j]);
					all_vals_key2.push_back(col_v_key2[j]);
				}

				MPI_Send(&length,1,MPI_INT,i,0,MPI_COMM_WORLD);
	 			MPI_Send(&col_v_countcum[i*size],size+1,MPI_INT,i,1,MPI_COMM_WORLD);
	 			MPI_Send(&all_vals_key1[0],length,MPI_INT,i,2,MPI_COMM_WORLD);	
	 			MPI_Send(&all_vals_key2[0],length,MPI_FLOAT,i,3,MPI_COMM_WORLD);

			}
		}

		vectlength=0;
		if (world_rank!=0){
		  	MPI_Recv(&vectlength,1, MPI_INT, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		vector < int > v4(size+1,0);
		vector < int > v5(vectlength,0);
		vector < float > v6(vectlength,0);
		if (world_rank!=0){
			MPI_Recv(&v4[0],size+1, MPI_INT, 0,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&v5[0],vectlength, MPI_INT, 0,2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  	MPI_Recv(&v6[0],vectlength, MPI_FLOAT, 0,3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			vector< std::pair <int,float> > vect;
			for (int i=0;i<vectlength;i++){
				vect.push_back(make_pair(v5[i],v6[i]));
			}
			
			int indexstart = v4[0];
			for (int i=0;i<size;i++){
				vector< std::pair <int,float> > keys;
				int start = v4[i];
				int end = v4[i+1];
				for (int j=start-indexstart;j<end-indexstart;j++){
					keys.push_back(vect[j]);
				}
				stable_sort(keys.begin(),keys.end(),[](const std::pair<int, float>& p1, const std::pair<int, float>& p2) { return p1.first < p2.first;});	
				for (int j=0;j<keys.size();j++){
					vect[start-indexstart+j] = keys[j];
				}	
			}
			for (int i=0;i<vect.size();i++){
				v5[i]=vect[i].first;
				v6[i]=vect[i].second;
			}
			MPI_Send(&vectlength,1,MPI_INT,0,0,MPI_COMM_WORLD);
			MPI_Send(&v5[0],vectlength,MPI_INT,0,1,MPI_COMM_WORLD);	
	 		MPI_Send(&v6[0],vectlength,MPI_FLOAT,0,2,MPI_COMM_WORLD);
		}


		if (world_rank==0){

			for (int i=1;i<world_size;i++){
				int vectlength;
		  		MPI_Recv(&vectlength,1, MPI_INT, i,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				vector<int> all_vals_key1(vectlength,0);
	 			vector<float> all_vals_key2(vectlength,0);
				MPI_Recv(&all_vals_key1[0],vectlength, MPI_INT, i,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  		MPI_Recv(&all_vals_key2[0],vectlength, MPI_FLOAT, i,2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				int startpoint= col_v_countcum[i*size];				
				for (int j=0;j<vectlength;j++){
					col_v_key1[startpoint+j]=all_vals_key1[j];
					col_v_key2[startpoint+j]=all_vals_key2[j];
				}		
			}
			
			
			// SORT ITS OWN PROCESS DATA
			//INITIAL ROWS
			for (int i=0;i<size;i++){
				int start = col_v_countcum[i];
				int end = col_v_countcum[(i+1)];
				int length = end-start;
				vector< std::pair <int,float> > vect;
				for (int j=start;j<end;j++){
					vect.push_back(make_pair(col_v_key1[j],col_v_key2[j]));
				}
				stable_sort(vect.begin(),vect.end(),[](const std::pair<int, float>& p1, const std::pair<int, float>& p2) { return p1.first < p2.first;});	
				for (int j=0;j<vect.size();j++){
					col_v_key1[start+j]=vect[j].first;
					col_v_key2[start+j]=vect[j].second;
				}
			}

			//END ROWS
			for (int i=world_size*size;i<maxcols;i++){
				int start = col_v_countcum[i];
				int end = col_v_countcum[(i+1)];
				int length = end-start;
				vector< std::pair <int,float> > vect;
				for (int j=start;j<end;j++){
					vect.push_back(make_pair(col_v_key1[j],col_v_key2[j]));
				}
				stable_sort(vect.begin(),vect.end(),[](const std::pair<int, float>& p1, const std::pair<int, float>& p2) { return p1.first < p2.first;});	
				for (int j=0;j<vect.size();j++){
					col_v_key1[start+j]=vect[j].first;
					col_v_key2[start+j]=vect[j].second;
				}
			}

			
			for (int i=0;i<maxrows;i++){
				for (int j=0;j<maxcols;j++){
					int start = col_v_countcum[j];
					int end = col_v_countcum[(j+1)];
					int length = end-start;
					for (int k=0;k<length;k++){
						if (col_v_index[start+k]==i){
							row_v_key1.push_back(col_v_key1[start+k]);
							row_v_key2.push_back(col_v_key2[start+k]);
							row_v_index.push_back(j);
						}
					}
				}
			}


			col_v_index.clear();
			col_v_key1.clear();
			col_v_key2.clear();
			col_v_countcum.clear();	
		}
	}


	// Now print to file by process 0 by fwriting col_v_key and col_v_value_swapped in pairs

	if (world_rank==0){
		FILE *fp2;	
	    	fp2 = fopen(outputfile.c_str(), "wb"); 

		if(fp2 == NULL)
	    	{
	        	printf("Error opening file\n");
	        	exit(1);
	    	}
		int count=0;
	   	for (int i=0;i<maxrows;i++){
			int start = row_v_countcum[i];
			int end = row_v_countcum[(i+1)];
			int length = end-start;
			for (int l=0;l<length;l++){
				int row = i;
				int col = row_v_index[start+l];
				int val1 = row_v_key1[start+l];
				float val2 = row_v_key2[start+l];
				fwrite(&row,sizeof(int),1,fp2);
				fwrite(&col,sizeof(int),1,fp2);
				fwrite(&val1,sizeof(int),1,fp2);
				fwrite(&val2,sizeof(float),1,fp2);
				//cout<<row<<","<<col<<","<<val1<<","<<val2<<endl;
				count+=1;
				//if (col==218){cout<<row<<","<<col<<","<<val1<<","<<val2<<endl;}
				//if (j==218){cout<<i<<","<<j<<","<<row_v_key1[i*maxcolumns+j]<<","<<row_v_key2[i*maxcolumns+j]<<endl;}
			}
			
		}
		fclose (fp2);

	 }


		MPI_Finalize();
}


