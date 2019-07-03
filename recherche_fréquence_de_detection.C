#include <cmath>
#include <iostream>
#include <TMath.h>
#include <cmath>
#include "TF2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TGraph.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

void lire()
{
	double boss;	
	double m;	
	int compte = 10;	
	double c;	
	double j;	
	double resume[compte];	
	double r[compte];	
	int y = 1;	
	TFile* hz = new TFile("fichier_h=0Hz","READ");
	TTreeReader gh("Caracteristique", hz);
       	TTreeReaderValue<double_t> STD_0hz(gh, "Std");	
	TTreeReaderValue<double_t> Mean_0hz(gh, "Mean");
	for(int t = 0 ; t<compte ; t++)
	{	
				
		gh.Next();	
		TString compteur = "fichier_r=";		
		compteur +=y;				
		TFile* fileTree = new TFile(compteur,"READ");
		TTreeReader jk("Caracteristique", fileTree);
		TTreeReaderValue<double_t> Mean(jk, "Mean");
		jk.Next();
		r[t] = double(y)/10;		
		y = y + 1;		
		for (int i = 0 ; i < 60; i++)
		{
			if (i == 0)
			{
				boss = (*Mean);
				jk.Next();
			}			
			if (i != 0)
			{
				m = (*STD_0hz+*Mean_0hz)*5 - (*Mean);
				c = (*STD_0hz+*Mean_0hz)*5 - boss;		
				if ( TMath::Abs(m) < TMath::Abs(c))
				{
					//cout<<"le nouveau plus petit est: " <<(*Mean)<<endl;					
					boss = (*Mean);
					jk.Next();
				        j = i;
				}				
			
				if ( TMath::Abs(m) > TMath::Abs(c))
				{	
					jk.Next();
				}

			}
		}
		TTreeReader jl("Caracteristique", fileTree);		
		TTreeReaderValue<double_t> frequence(jl, "FrequenceSource");
		jl.Next();
		for (int p = 0 ; p < j ; p++)
		{		
			jl.Next();	
		
		}		
		cout<<"la frequence de detection est de: "<<*frequence<<"hz"<<" pour r="<<r[t]<<"°"<<endl;
		resume[t] = *frequence;
	}
        TCanvas* canvascercle = new TCanvas("frequence_detectable_en_fonction_de_la_taille_de_la_source","frequence_detectable_en_fonction_de_la_taille_de_la_source");
	TGraph* gr = new TGraph(compte,r,resume);       
	gr -> SetTitle("frequence de detection en focntion de la taille de la source"); 	
	gr->Draw("A*");	
	
}
