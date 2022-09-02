#include "Msh_reader.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstring> 
#include <cstdlib>
#include <sstream>

#include <cmath>
using namespace std;

void msh_read(std::string filename,double ** &vnodes,int * &vflags,int ** &velems,int &no_nodes,int &no_elements){
   
    ifstream file(filename.c_str());    
   
    
    //wczytanie pierwszego stringa
    char line[LINE_LENGTH];
    file.getline(line,LINE_LENGTH);
    
    //wczytanie liczby wezlow
    file.getline(line,LINE_LENGTH);
    stringstream iline;
    iline << line;
    iline.clear();
    iline >> no_nodes;
    iline.clear();

    if(msh_debug) {
    cout << "Liczba wezlow: " << no_nodes << endl;
    }

    //wczytanie wierzcholkow i flag
    double **nodes;
    int *flags;
    int id;
    nodes = new double*[no_nodes];
    for(int i=0;i<no_nodes;i++){
        nodes[i]=new double[3];
    }
    flags = new int[no_nodes];
    for(int i=0;i<no_nodes;i++){
        stringstream sline;
        file.getline(line,LINE_LENGTH);
        sline << line;
        sline >> id >> nodes[i][0] >> nodes[i][1] >> nodes[i][2];
    }
    
    
    file.getline(line,LINE_LENGTH);
    file.getline(line,LINE_LENGTH);
    
    //wczytyje liczbe elementow
    file.getline(line,LINE_LENGTH);
    iline << line;
    //int no_elements;
    iline >> no_elements;
    if(msh_debug) {
     cout << "Laczna liczba elementow fizycznych: " << no_elements << endl;
    }
     int **elements;
     elements = new int*[no_elements];
     for(int i=0;i<no_elements;i++){
         elements[i]=new int[4];
     }
     
    for(int i=0;i<no_nodes;i++){
	flags[i]=flag_interior;
    }

    //wczytuje elementy
     int ele_iter = 0;
     for(int i=0;i<no_elements;i++){
        stringstream sline;
        int tmp_int[9];
        file.getline(line,LINE_LENGTH);
        sline << line;
        // w zaleznosci od rodzaju flagi wykonywane beda rozne czynnosci
        sline >> tmp_int[0] >> tmp_int[1] >> tmp_int[2] >> tmp_int[3] >> tmp_int[4] ;
        
        switch( tmp_int[1] ){
            case flag_point: break;
            case flag_line: //break;
			switch( tmp_int[2] ){
				case flag_d1:
				case flag_d2: 
				case flag_d3:
				case flag_d4:
				case flag_d5:
				case flag_n1: 
					if( tmp_int[4] == 1 ){
						sline >> tmp_int[5] ;
						flags[ tmp_int[5]-1 ] = tmp_int[2];
						   if(msh_debug) cout << "F[" << i << "]=" << tmp_int[5] << endl;
					}else if (tmp_int[4] == 2) {
						   sline >> tmp_int[5]  >> tmp_int[6];
						   flags[ tmp_int[5]-1 ] = tmp_int[2];
						   flags[ tmp_int[6]-1 ] = tmp_int[2];
						   if(msh_debug) cout << "F[" << i << "]=" << tmp_int[5] <<", " << tmp_int[6] << endl;
					}else if (tmp_int[4] == 3){
						   sline >> tmp_int[5]  >> tmp_int[7];
						   if(msh_debug) cout << "F[" << i << "]=" << tmp_int[5] <<", " << tmp_int[6] << endl;
					}else{
						   cout << "<Blad> Czy plik z biliniowymi elementami?" << endl;
					}
					break;
			}
			break;
            //~ case flag_interior:
            case flag_triangle:
                sline >> tmp_int[5] >> tmp_int[6] >> tmp_int[7];// >> tmp_int[8] ; // zakomentowane bo to sa trojkaty teraz
                for(int j=0;j<4;j++){
                        elements[ele_iter][j]=tmp_int[j+5]; 
                }
                if(msh_debug) cout <<  "E[" << ele_iter << "]= " << elements[ele_iter][0] << ", " << elements[ele_iter][1] << ", " << elements[ele_iter][2] << endl;
                ele_iter = ele_iter + 1;

                break;
            default:
                cout << "<Blad> Flaga nie znana: " << tmp_int[2] << endl;
        }
        
    }
     
     //alokujemy nasze tablice
    vnodes = new double*[no_nodes];
    for(int i=0;i<no_nodes;i++){
        vnodes[i]=new double[3];
    }
    
     velems = new int*[ele_iter];
     for(int i=0;i<ele_iter;i++){
         velems[i]=new int[4];
     }


    vflags = new int[no_nodes];
    
    //kopiowanie do naszych tablic
    for(int i=0;i<no_nodes;i++){
        for(int j=0;j<3;j++){
            vnodes[i][j]=nodes[i][j];
        }
        vflags[i]=flags[i];
    }
    
    for(int i=0;i<ele_iter;i++){
        for(int j=0;j<4;j++){
            velems[i][j]=elements[i][j]-1;
        }
    }

    if(msh_debug) {
        cout <<"Wezly:" <<endl;
        for(int i=0;i<no_nodes;i++){
            cout << i <<": " << flags[i] << ", " << nodes[i][0] << ", " << nodes[i][1] << ", " << nodes[i][2] << ", " << flags[i] << endl;
        }
    }

    //dealokacja
    for(int i=0;i<no_nodes;i++){
        delete [] nodes[i];
    }
    delete [] nodes;
    delete [] flags;
    for(int i=0;i<no_elements;i++){
        delete [] elements[i];
    }
    delete [] elements;
    
    
    file.close();  

    no_elements=ele_iter;
}


void deallocate_tables(double **vnodes,int *vflags,int no_nodes){
    
    for(int i=0;i<no_nodes;i++){
        delete [] vnodes[i];
    }
    delete [] vnodes;
    delete [] vflags;
    
}


