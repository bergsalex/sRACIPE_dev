#include <Rcpp.h>
#include <iostream>
#include <random>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>

using namespace Rcpp;


// [[Rcpp::export]]

IntegerMatrix interaction_reader(IntegerMatrix gene_interaction, String filepath, String filename, int number_gene)
  {

  std::string file_name = filename;
  std::fstream out_2("./results/gene_interaction_topology_" + file_name + ".txt",std::fstream::out);
  file_name = filepath;
  //file_name="./inputs/"+file_name;
  std::ifstream infile ( file_name, std::ifstream::in);
  if(!infile) {
    Rcout << "Cannot open input file.\n";
    return 1;
  }
  std::vector <std::string> gene_names;
  std::string word;
  int i=0; int first_gene=0,second_gene=0;
  while (infile >> word)
  {
    if(i<3){}
    else if(i%3!=2)
    {
      if (std::find(gene_names.begin(), gene_names.end(),word)!=gene_names.end())
      {
        if(i%3==0)
        {first_gene = distance(gene_names.begin(), find(gene_names.begin(), gene_names.end(), word));}
        else second_gene = distance(gene_names.begin(), find(gene_names.begin(), gene_names.end(), word));
      }
      else
      {
        gene_names.push_back(word);
        if((i%3==0))first_gene=gene_names.size()-1;
        else second_gene=gene_names.size()-1;
      }

    }
    else gene_interaction(second_gene,first_gene)=std::stoi(word);
    i++;
  }

  out_2<<"gene"<<"\n";
  for (int i=0; i<gene_names.size(); i++)
  {
    out_2 <<gene_names[i] << "\n";
  }
      /*
  out_2<<"\n";
    for (int i=0; i<gene_names.size(); i++)
  {
        out_2 <<gene_names[i] << "\t";
      for (int j=0; j<gene_names.size(); j++)
        out_2<<gene_interaction(i,j)<<"\t";
    out_2<<"\n";

  }
       */
  Rcout<<gene_names.size()<<"\n";

  return gene_interaction;

}
