#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libxml/xmlwriter.h>
#include <libxml/xmlstring.h>

#include "allvars.h"
#include "proto.h"

struct io_header_1 header1;

int NumPart, Ngas;

struct particle_data
{
  double Pos[3];
  double Vel[3];
  double Mass;
  int Type;

  double Rho, U, Temp, Ne;
} *Pd;

int *Id;

double Time, Redshift;



void rw_snapshot(void){
  char input_fname[200];


  sprintf(input_fname, "%s/%s", OutputDir, FileBase);
  printf("Reading %d Files \n",NumFilesWrittenInParallel+1);

  load_snapshot(input_fname, NumFilesWrittenInParallel+1);


  reordering();			/* call this routine only if your ID's are set properly */

  unit_conversion();		/* optional stuff */

  write_data(input_fname);


}





/* here we write the xml file
 */
void write_data(char *fname)
{
	char buf[200];

	sprintf(buf, "%s.xml", fname);
	xmlTextWriterPtr writer = xmlNewTextWriterFilename(buf, 0);
	            xmlTextWriterStartDocument(writer, NULL, NULL, NULL);
	            {
	                xmlTextWriterStartElement(writer, BAD_CAST "N-GenIC");
	                {
					int i;
			//		sprintf(buf, "%s.txt", fname);
					//FILE *fd = fopen (buf,"w");
					printf("\n Saving xml...");
					 sprintf(buf, "%d", NumPart);
				 	 xmlTextWriterWriteAttribute(writer, BAD_CAST "NumBodies", xmlCharStrdup(buf));
					 sprintf(buf, "%f", Redshift);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Redshift", xmlCharStrdup(buf));
					 sprintf(buf, "%f", Box/1000);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Box", xmlCharStrdup(buf));
					 sprintf(buf, "%f", Omega);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "OmegaMatter", xmlCharStrdup(buf));
					 sprintf(buf, "%f", OmegaLambda);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "OmegaLambda", xmlCharStrdup(buf));
					 sprintf(buf, "%f", OmegaBaryon);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "OmegaBaryon", xmlCharStrdup(buf));
					 sprintf(buf, "%f", 10.0);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "HubbleParam_10km_s_Mpc_h", xmlCharStrdup(buf));
					 if(WhichSpectrum==1)
						 sprintf(buf, "Eisenstein and Hu spectrum");
					 else if(WhichSpectrum==0)
						 sprintf(buf, "Efstathiou parametrization");
					 else
						 sprintf(buf, "Tabulated power spectrum");

					 xmlTextWriterWriteAttribute(writer, BAD_CAST "WhichSpectrum", xmlCharStrdup(buf));
					 sprintf(buf, "%f", Sigma8);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Sigma8", xmlCharStrdup(buf));
					 sprintf(buf, "%f", ShapeGamma);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "ShapeGamma", xmlCharStrdup(buf));
					 sprintf(buf, "%f",PrimordialIndex);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "PrimordialIndex", xmlCharStrdup(buf));
					 sprintf(buf, "%e", UnitLength_in_cm/100*1000/HubbleParam);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "UnitLength_in_m_Mpc_h", xmlCharStrdup(buf));
					 sprintf(buf, "%e", UnitMass_in_g/1000/HubbleParam);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "UnitMass_in_kg_1e10_Msol_h", xmlCharStrdup(buf));
					 sprintf(buf, "%e",  10000.0);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "UnitVelocity_in_m_s_10Km_s", xmlCharStrdup(buf));
					 sprintf(buf, "%e", GRAVITY/1000);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Gravitational_Constant_m3_kg_s2", xmlCharStrdup(buf));
					 sprintf(buf, "Yes");
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Binary_file", xmlCharStrdup(buf));
					 sprintf(buf, "%s.bin", fname);
					 xmlTextWriterWriteAttribute(writer, BAD_CAST "Data_filename", xmlCharStrdup(buf));

					 printf("...done");
					 printf("\nSaving ICs...");

				   FILE* particles = fopen(buf, "wb");
				   FILE* fd = fopen("ics.txt", "w");

				   for(i = 1; i <= NumPart; i++)
						{

						  Pd[i].Pos[0]=Pd[i].Pos[0]/1000; Pd[i].Pos[1]=Pd[i].Pos[1]/1000; Pd[i].Pos[2]=Pd[i].Pos[2]/1000;
					  Pd[i].Vel[0]=Pd[i].Vel[0]/10; Pd[i].Vel[1]=Pd[i].Vel[1]/10; Pd[i].Vel[2]=Pd[i].Vel[2]/10;
				 fprintf(fd,"%4.16f %4.16f %4.16f %4.16f %4.16f %4.16f %4.16f\n",Pd[i].Mass,Pd[i].Pos[0],Pd[i].Pos[1],Pd[i].Pos[2],Pd[i].Vel[0],Pd[i].Vel[1],Pd[i].Vel[2]);
					    fwrite ((void*)&Pd[i].Mass, 1 , sizeof(double) , particles );
					    fwrite((void*)&Pd[i].Pos[0], 3*sizeof(double), 1, particles);
					    fwrite((void*)&Pd[i].Vel[0], 3*sizeof(double), 1, particles);
				// printf("\n%f %f %f %f %f %f %f\n",Pd[i].Mass,Pd[i].Pos[0],Pd[i].Pos[1],Pd[i].Pos[2],Pd[i].Vel[0],Pd[i].Vel[1],Pd[i].Vel[2]);
						}

					printf("...done.\n");
					fclose(fd);

					xmlTextWriterEndElement(writer); // particles
				   fclose(particles);
	                }
				xmlTextWriterEndElement(writer);
	            }
			xmlTextWriterEndDocument(writer);
			xmlFreeTextWriter(writer);
   Pd++;
  free(Pd);

}





/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
void unit_conversion(void)
{
  //double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;

  /* physical constants in cgs units */
  //GRAVITY = 6.672e-8;
  //BOLTZMANN = 1.3806e-16;
  //PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm = 3.0856776e24;	/*  code length unit in cm */
  UnitMass_in_g = 1.9885e43;	/*  code mass unit in g */
  UnitVelocity_in_cm_per_s = 1.0e5;

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);


  Xh = 0.76;			/* mass fraction of hydrogen */
  HubbleParam = 0.65;


  for(i = 1; i <= NumPart; i++)
    {
      if(Pd[i].Type == 0)	/* gas particle */
	{
	  MeanWeight = 4.0 / (3 * Xh + 1 + 4 * Xh * Pd[i].Ne) * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u = Pd[i].U * UnitEnergy_in_cgs / UnitMass_in_g;

	  gamma = 5.0 / 3;

	  /* get temperature in Kelvin */

	  Pd[i].Temp = MeanWeight / BOLTZMANN * (gamma - 1) * u;
	}
    }

}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
void load_snapshot(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i, k, dummy, ntot_withmasses;
  int n, pc, pc_new, pc_sph,flag;

#define SKIP flag=fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
       sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      printf("reading heater `%s' ...\n", buf);
      fflush(stdout);

      flag=fread(&dummy, sizeof(dummy), 1, fd);
      flag=fread(&header1, sizeof(header1), 1, fd);
      flag=fread(&dummy, sizeof(dummy), 1, fd);
      if(files == 1)
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    NumPart += header1.npart[k];
	  Ngas = header1.npart[0];
	}
      else
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
		{NumPart += header1.npartTotal[k];
	  Ngas = header1.npartTotal[0];}
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1.mass[k] == 0)
	    ntot_withmasses += header1.npart[k];
	}

      if(i == 0)
	allocate_memory();

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
		  flag=fread(&Pd[pc_new].Pos[0], sizeof(double), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
		  flag=fread(&Pd[pc_new].Vel[0], sizeof(double), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;


      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
		  flag=fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses > 0)
	SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      Pd[pc_new].Type = k;

	      if(header1.mass[k] == 0)
	    	  flag=fread(&Pd[pc_new].Mass, sizeof(double), 1, fd);
	      else
		Pd[pc_new].Mass = header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses > 0)
	SKIP;


      if(header1.npart[0] > 0)
	{
	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
		  flag=fread(&Pd[pc_sph].U, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
		  flag=fread(&Pd[pc_sph].Rho, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
		{
	    	  flag=fread(&Pd[pc_sph].Ne, sizeof(double), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	      {
		Pd[pc_sph].Ne = 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }


  Time = header1.time;
  Redshift = header1.redshift;
}




/* this routine allocates the memory for the 
 * particle data.
 */
void allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(Pd = malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  Pd--;				/* start with offset 1 */


  if(!(Id = malloc(NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  Id--;				/* start with offset 1 */

  printf("allocating memory...done\n");
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
void reordering(void)
{
  int i;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i = 1; i <= NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource = Pd[i];
	  idsource = Id[i];
	  dest = Id[i];

	  do
	    {
	      psave = Pd[dest];
	      idsave = Id[dest];

	      Pd[dest] = psource;
	      Id[dest] = idsource;

	      if(dest == i)
		break;

	      psource = psave;
	      idsource = idsave;

	      dest = idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;
  free(Id);

  printf("space for particle ID freed\n");
}
