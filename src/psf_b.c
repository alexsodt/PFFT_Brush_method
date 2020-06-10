#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "psf_b.h"

/**
 * Read the total atom count from a PSF.
 *
 * Retrieves the total number of atoms from the NATOM field of a PSF. This
 * method does not verify that the NATOM field matches the number of atoms
 * listed. If the file is unreadable, or the NATOM field cannot be read, a value
 * of -1 is returned instead.
 *
 * @param[in] path The PSF from which the atom count should be read.
 * @return The number of atoms in the PSF.
 */
int getNATOM(char * path) {
  FILE * psf = fopen(path,"r");
  if(!psf) // Error encountered while opening file.
    return -1;
  int natom = -1; // in case no atoms are read successfully.
  char bufferA[1024]; // previous word
  char bufferB[1024]; // current word
  int numread = -1;
  while ( (!ferror(psf)) && (!feof(psf)) && numread!=0) { //stop looping upon error/eof
    if(!strncmp(bufferB,"!NATOM",5)) { //watch for !NATOM
      natom = atoi(bufferA); // read word (number) before !NATOM
      break; //stop looping
    }
    strcpy(bufferA,bufferB); //save current word as previous
    numread = fscanf(psf," %1024s ",bufferB); //load next word
  }
  fclose(psf);
  return natom;
}

/**
 * Read scattering lengths from a PSF.
 *
 * Reads scattering lengths from a PSF file and stores the contents in the
 * provided array. The format should be a standard PSF with an appended
 * scattering length column  (SLB), with the appropriate first-line signature:
 * PSF SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,1X,G14.6)
 * PSF XPLOR SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,1X,G14.6)
 * PSF CMAP CHEQ SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6,1X,G14.6)
 * PSF CMAP CHEQ XPLOR SLB
 * (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6,1X,G14.6)
 * PSF EXT SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,1X,G14.6)
 * PSF EXT XPLOR SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,1X,G14.6)
 * PSF EXT CMAP CHEQ SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6,1X,G14.6)
 * PSF EXT CMAP CHEQ XPLOR SLB
 * (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6,1X,G14.6)
 * The general format is:
 * II,LSEGID,LRESID,LRES,TYPE,IAC,CG,AMASS,IMOVE,[ECH,EHA],SLD
 * Formats sourced from charmm src file: source/io/psfres.src
 * The DRUDE extension to the PSF format is NOT supported, simply because the
 * impact of the DRUDE extension on a PSF is not understood by this author.
 *
 * @param[in] path The PSF from which the scattering lengths will be read.
 * @param[out] bs The array into which the scattering lengths will be stored.
 * @return The number of scattering lengths read, or zero if an error occurs
 * before any values are read.
 */
int psfToBs(char * path, double * bs) {
  FILE * psf = fopen(path,"r");
  if(!psf) // Error encountered while opening file.
    return 0;

  char buffer[1024]; //For storing one line of the file at a time.

  int offset=71; // Base offset for PSF

  char * ptr;

  ptr = fgets(buffer,1024,psf); //get first line
  if(ferror(psf) || feof(psf) || ptr==NULL) {
    fclose(psf);
    return 0; // Error reading first line of file.
  }
  if(strncmp(buffer,"PSF",3)) {
    fclose(psf);
    return 0; // File is not a PSF.
  }

  // Processing first line.
  if(!strstr(buffer,"SLB")) {
    fclose(psf);
    return 0; // File does not contain scattering lengths.
  }
  if(strstr(buffer,"EXT")) {
    offset+=18; // Adjust offset appropriately for extended PSF.
    if(strstr(buffer,"XPLOR"))
      offset+=2; // expanded XPLOR atom types get an extra two characters.
  }
  if(strstr(buffer,"CMAP CHEQ"))
    offset+=28; // Adjust offset to account for the extra two columns.

  //Loop until buffer contains the !NATOM line.
  while (1) {
    ptr = fgets(buffer,1024,psf); //get next line
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      fclose(psf);
      return 0; //Error reading file.
    }
    if(strstr(buffer,"!NATOM")) {
      break; //Found the atoms section.
    }
  }

  //Start reading scattering lengths
  char slbstring[15]; // For storing the scattering length as a string.
  slbstring[14]='\0'; // We'll only copy 14 characters, so the 15th needs to be
                      // the null terminator.
  int n = 0; //track the number of scattering lengths read so far.
  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      fclose(psf);
      return n; // Error reading file for next scattering length.
    }
    if(buffer[0]=='\n') {
      fclose(psf);
      return n; // Reached the end of the atoms section.
    }
    if(strstr(buffer,"!")) {
      fclose(psf);
      return n; // PSF is missing empty line between sections.
    }
    memcpy(slbstring,&buffer[offset],14); //isolate the SLB substring.
    bs[n]=atof(slbstring); //add the scatering length to the provided array.
    n++;
  }
}

/**
 *
 * Read residue names from a PSF.
 *
 * Reads residue names from a PSF file and stores the contents in the
 * provided array.
 *
 * @param[in] path The PSF from which the scattering lengths will be read.
 * @param[out] resnames The array into which the residue names will be stored.
 * @return The number of residue names read, or zero if an error occurs
 * before any values are read.
 */
int psfToResNames(char * path, char resnames[][9]) {
  FILE * psf = fopen(path,"r");
  if(!psf) // Error encountered while opening file.
    return 0;

  char buffer[1024]; //For storing one line of the file at a time.

  int offset=19; // Base offset for PSF
  int resnamewidth=4; //Default for non-extended PSF

  char * ptr;

  ptr = fgets(buffer,1024,psf); //get first line
  if(ferror(psf) || feof(psf) || ptr==NULL) {
    fclose(psf);
    return 0; // Error reading first line of file.
  }
  if(strncmp(buffer,"PSF",3)) {
    fclose(psf);
    return 0; // File is not a PSF.
  }

  // Processing first line.
  if(strstr(buffer,"EXT")) {
    offset+=10; // Adjust offset appropriately for extended PSF.
    resnamewidth=8; //extedned PSF files have 8-character residue names.
  }

  //Loop until buffer contains the !NATOM line.
  while (1) {
    ptr = fgets(buffer,1024,psf); //get next line
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      fclose(psf);
      return 0; //Error reading file.
    }
    if(strstr(buffer,"!NATOM")) {
      break; //Found the atoms section.
    }
  }

  //Start reading scattering lengths
  char resname[resnamewidth+1]; // For storing the scattering length as a string.
  resname[resnamewidth]='\0'; // We'll only copy 4 or 8 characters, so the 5th or 9th needs to be
                      // the null terminator.
  int n = 0; //track the number of residue names read so far.
  while (1) {
    ptr = fgets(buffer,1024,psf);
    if(ferror(psf) || feof(psf) || ptr==NULL) {
      fclose(psf);
      return n; // Error reading file for next scattering length.
    }
    if(buffer[0]=='\n') {
      fclose(psf);
      return n; // Reached the end of the atoms section.
    }
    if(strstr(buffer,"!")) {
      fclose(psf);
      return n; // PSF is missing empty line between sections.
    }
    memcpy(resname,&buffer[offset],resnamewidth); //isolate the SLB substring.
    strcpy(resnames[n],resname); //add the scatering length to the provided array.
    n++;
  }
}

/**
 * Just a quick main method for testing. Takes a psf file as its only argument.
 */
//int main( int argc, char **argv) {
//  int natom = getNATOM(argv[1]);
//  printf("%d\n",natom);
//  double bs[natom];
//  int n = psfToBs(argv[1], bs);
//  int i;
//  
//  //printing 10 per line
//  int wholelines = natom/10;
//  int remainder = natom % 10;
//  if (remainder==0) { // the last line is special, so make sure we have one.
//    wholelines--;
//    remainder+=10;
//  }
//  for (i=0; i<wholelines; i++) {
//    for (int j=0; j<10; j++)
//      printf("%lf,",bs[10*i+j]);
//    printf("\n");
//  }
//  for (int j=0; j<(remainder-1); j++)
//    printf("%lf,",bs[10*wholelines+j]);
//  printf("%lf\n",bs[10*wholelines+remainder-1]);
//  if((natom-1) != 10*wholelines+remainder-1)
//    printf("SOMETHING WENT WRONG!\n");
//}
