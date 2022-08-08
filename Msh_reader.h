
#ifndef MSH_READER_H
#define	MSH_READER_H
#include <iostream>
#include <string>
#define LINE_LENGTH 200

#define msh_debug 0
#define flag_d1  11 // dirichlety
#define flag_d2  22 
#define flag_d3  33  
#define flag_d4  44 
#define flag_d5  55  
#define flag_n1 166 // neumany
#define flag_interior   13
#define flag_point   15
#define flag_line   1
#define flag_triangle   2
        
    /*
     * Wczytuje z pliku o nazwie filename wspolrzedne wezlow, flagi oraz elementy, a także liczbe wezlow.
     * @param filename nazwa pliku z danymi
     * @param vnodes talica, do ktorej zapisane sa wspolrzedne x,y,z wezlow. Wymiar tablicy po wyjsciu z funkcji [no_nodes][3]
     * @param vflags tablica, do której zapisane sa flagi wezlow 
     * @param no_nodes liczba wczytanych wezlow
     * 
     */
    //void msh_read(std::string filename,double ** &vnodes,int * &vflags,int & no_nodes);
    void msh_read(std::string filename,double ** &vnodes,int * &vflags,int ** &velements,int &no_nodes,int &no_elements);
    /*
     * Dealokuje tablice zaalokowane w readFile. 
     * @param vnodes talica, do ktorej zapisane sa wspolrzedne x,y,z wezlow. Wymiar tablicy po wyjsciu z funkcji [no_elements][3]
     * @param vflags tablica, do której zapisane sa flagi wezlow 
     * @param no_nodes liczba wczytanych wezlow
     * 
     */
    void deallocate_tables(double **vnodes,int *vflags,int no_nodes);


    
#endif	/* MSH_READER_H */
