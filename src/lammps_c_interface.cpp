//
// Created by olwhi on 24/07/2023.
//

#include "lammps_c_interface.h"
#include <stdio.h>
#include <cstdint>
#include <unistd.h>
#include <fstream>

LammpsObject::LammpsObject(){

}

LammpsObject::LammpsObject(int selector, string prefixIn, string prefixOut) {
    cout << "Creating Lammps Object" << endl;
    const char *lmpargv[] = { "liblammps", "-screen", "none"};
//    const char *lmpargv[] = { "liblammps", "-log", "none"};
    prefixFolderIn=prefixIn;
    prefixFolderOut=prefixOut;
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);

    cout << "lmpargc -> " << lmpargc << endl;
    
    /* create LAMMPS instance */
    handle = lammps_open_no_mpi(lmpargc, (char **)lmpargv, NULL);

    cout << "created lammps handle" << endl;

    if (handle == NULL) {
        printf("LAMMPS initialization failed");
        lammps_mpi_finalize();
        exit(1);
    }

    /* get and print numerical version code */
    version = lammps_version(handle);
    printf("LAMMPS Version: %d\n",version);
    string fname;
    if (selector==0){          cout << "-------------------------> Simple Graphene" << endl; fname=prefixIn+"/Si.in";      }
    else if (selector==1){     cout << "-------------------------> Triangle Raft  " << endl; fname=prefixIn+"/Si2O3.in";   }
    else if (selector==2){     cout << "-------------------------> Bilayer        " << endl; fname=prefixIn+"/SiO2.in";    }
    else if (selector==3){     cout << "-------------------------> Tersoff        " << endl; fname=prefixIn+"/C.in"; }
    else if (selector==4){     cout << "-------------------------> BN             " << endl; fname=prefixIn+"/BN.in"; }
    else{exit(3);}
    cout << "File to run from : " << fname << endl;
    lammps_file(handle, fname.c_str());

    natoms = (int)(lammps_get_natoms(handle)+0.5);
    nbonds = std::intptr_t (lammps_extract_global(handle, "nbonds"));
    nangles= std::intptr_t (lammps_extract_global(handle, "nangles"));

}




int LammpsObject::initialiseLammpsObject() {
    const char *lmpargv[] = { "liblammps", "-log", "none"};
    int lmpargc = sizeof(lmpargv)/sizeof(const char *);

    /* create LAMMPS instance */
    handle = lammps_open_no_mpi(lmpargc, (char **)lmpargv, NULL);
    if (handle == NULL) {
        printf("LAMMPS initialization failed");
        lammps_mpi_finalize();
        return 1;
    }

    /* get and print numerical version code */
    version = lammps_version(handle);
    printf("LAMMPS Version: %d\n",version);
    return 0;
}

int LammpsObject::write_data(int selector){
    string fname;
    if (selector==0){          fname="write_data "+prefixFolderOut+"/Si_results.dat";      }
    else if (selector==1){     fname="write_data "+prefixFolderOut+"/Si2O3_results.dat";   }
    else if (selector==2){     fname="write_data "+prefixFolderOut+"/SiO2_results.dat";    }
    else if (selector==3){     fname="write_data "+prefixFolderOut+"/C_results.dat";       }
    else if (selector==4){     fname="write_data "+prefixFolderOut+"/BN_results.dat";      }
    else{exit(3);}

    lammps_command(handle, fname.c_str());
    return 0;
}

int LammpsObject::write_restart(int selector){
    string fname;
    if (selector==0){          fname="write_restart "+prefixFolderOut+"/Si_restart.restart";      }
    else if (selector==1){     fname="write_restart "+prefixFolderOut+"/Si2O3_restart.restart";   }
    else if (selector==2){     fname="write_restart "+prefixFolderOut+"/SiO2_restart.restart";    }
    else if (selector==3){     fname="write_restart "+prefixFolderOut+"/C_restart.restart";       }
    else if (selector==4){     fname="write_restart "+prefixFolderOut+"/BN_restart.restart";      }
    else{exit(3);}

    lammps_command(handle, fname.c_str());

    return 0;
}

int LammpsObject::finaliseLammpsObject(int selector){
    string fname;
    if (selector==0){          fname="write_data "+prefixFolderOut+"/Si_results.dat";      }
    else if (selector==1){     fname="write_data "+prefixFolderOut+"/Si2O3_results.dat";   }
    else if (selector==2){     fname="write_data "+prefixFolderOut+"/SiO2_results.dat";    }
    else if (selector==3){     fname="write_data "+prefixFolderOut+"/C_results.dat";    }
    else if (selector==4){     fname="write_data "+prefixFolderOut+"/BN_results.dat";      }

    else{exit(3);}

    lammps_command(handle, fname.c_str());
    lammps_close(handle);
    lammps_mpi_finalize();
    return 0;
}

void LammpsObject::runInput(string fname){
    lammps_file(handle, fname.c_str());

//    natoms = std::intptr_t (lammps_extract_global(handle, "natoms"));
    nbonds = std::intptr_t (lammps_extract_global(handle, "nbonds"));
    nangles= std::intptr_t (lammps_extract_global(handle, "nangles"));

}

void LammpsObject::getatominfo(int dim){

    //std::intptr_t natom, nbond, nangle;


//    natoms = std::intptr_t (lammps_extract_global(handle, "natoms"));
    nbonds = std::intptr_t (lammps_extract_global(handle, "nbonds"));
    nangles= std::intptr_t (lammps_extract_global(handle, "nangles"));

    // gather atoms
//    x = new double[dim*natoms];
//    lammps_gather_atoms(handle, "x", 1, dim, x);

//    b = new double[3*nbonds];
//    lammps_gather_bonds(handle, b);

//    angles = new double[4*nangles];
//    lammps_gather_angles(handle, angles);

}
/*
double LammpsObject::fetchBonds(){
    nbonds = std::intptr_t (lammps_extract_global(handle, "nbonds"));
    bonds = new double[3*nbonds];
    lammps_gather_bonds(handle, bonds);
    return bonds;
}
*/
// actually unmanagable as order will change (?)
/*
 * bool LammpsObject::compareBonds(double* old_bonds){
    nbonds = std::intptr_t (lammps_extract_global(handle, "nbonds"));
    bonds = new double[3*nbonds];
    lammps_gather_bonds(handle, bonds);

}
*/


double LammpsObject::pbx(){
    return *(double *)lammps_extract_global(handle, "boxxhi");
}
double LammpsObject::pby(){
    return *(double *)lammps_extract_global(handle, "boxyhi");
}
double LammpsObject::pbz(){
    return *(double *)lammps_extract_global(handle, "boxzhi");
}

int LammpsObject::getnAtoms(){
    return natoms;
}

double* LammpsObject::fetchCrds(int dim){
//    natoms = std::intptr_t (lammps_extract_global(handle, "natoms"));
    // gather atoms
    coords= new double[dim*natoms];
    lammps_gather_atoms(handle, "x", 1, dim, coords);
    return coords;
}

void LammpsObject::pushCrds(int dim, double* old_coords){
    //natoms = std::intptr_t (lammps_extract_global(handle, "natoms"));
    // gather atoms
    //coords= new double[dim*natoms];
    lammps_scatter_atoms(handle, "x", 1, dim, old_coords);
}

double* LammpsObject::fetchBonds(){
    double initialBondCount= lammps_get_thermo(handle, "bonds");
    int bond_size = 3*initialBondCount;
    cout << bond_size << endl;
    bonds = new double[bond_size];
    lammps_gather_bonds(handle, bonds);
    return bonds;
}

void LammpsObject::breakBond(int atom1, int atom2, int type) {
    double initialBondCount= lammps_get_thermo(handle, "bonds");

    string command;
    if (globalVerbose){cout << "                                                              delete bond " << atom1 <<" "<< atom2 << endl;}
    command = "group switch id " + to_string(atom1) + " " + to_string(atom2);
    lammps_command(handle, command.c_str());
    command = "delete_bonds switch bond " + to_string(type) + " remove";
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");
//    lammps_command(handle, "run 0");
    double finalBondCount = lammps_get_thermo(handle, "bonds");
    if (globalVerbose){cout << "Initial " << initialBondCount << " vs Final " << finalBondCount << endl;}
    if (finalBondCount!=initialBondCount-1){
        cout << "Error in Bond Counts" << endl;
        exit(7);
    }

}

void LammpsObject::formBond(int atom1, int atom2, int type){
    string command;
    if (globalVerbose){cout << "                                                              create bond " << atom1 <<" "<< atom2 << endl;}
    command = "create_bonds single/bond " +to_string(type) + " " + to_string(atom1) + " " + to_string(atom2);
    lammps_command(handle, command.c_str());
}

void LammpsObject::breakAngle(int atom1, int atom2, int atom3){
    double initialAngleCount= lammps_get_thermo(handle, "angles");


    string command;
    //cout << globalVerbose << endl;
    if (globalVerbose) {cout << "                                                              delete angle " << atom1 <<" "<< atom2 << " "<< atom3 << endl;}

    command = "group switch id " + to_string(atom1) + " " + to_string(atom2) + " " + to_string(atom3);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");
//    lammps_command(handle, "run 0");
    double finalAngleCount = lammps_get_thermo(handle, "angles");
    if (globalVerbose) {cout << "Initial " << initialAngleCount << " vs Final " << finalAngleCount << endl;}
    if (finalAngleCount!=initialAngleCount-1){
        cout << "Error in Angle Counts" << endl;
        exit(7);
    }

}

void LammpsObject::formAngle(int atom1, int atom2, int atom3) {
    string command;
    if (globalVerbose) {cout << "                                                              create angle " << atom1 <<" "<< atom2 << " "<< atom3 << endl;}
    command = "create_bonds single/angle 1 " + to_string(atom1) + " " + to_string(atom2) + " " + to_string(atom3);
    lammps_command(handle, command.c_str());
}


// solution here likely involves grabbing angles and bonds from Network
void LammpsObject::switchGraphene(VecF<int> switchIdsA, Network networkA) {
    if (globalVerbose) cout << "----------------- Switch Graphene" << endl;
    bool verbose=false;
    int nbonds0, nbonds1, nangles0, nangles1;

    //unpck parameters
    int a,b,c,d,e,f;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];


    a +=1;
    b +=1;
    c +=1;
    d +=1;
    e +=1;
    f +=1;
    /*
    if (a==266 || b==266 || c==266 || d==266 || e==266 || f==266 || a==227 || b==227 || c==227 || d==227 || e==227 || f==227) {

        cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;

//        307 267 347 308 266 228
//        266 227 306 186 267 187

//        266 227 306 186 267 187


        sleep(10);
        cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n INCORPORATED 267 \n\n\n" << endl;
    }
    */
    if (verbose) {


        cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;
        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[a - 1].netCnxs[i] + 1 << " ";
        }
        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[b - 1].netCnxs[i] + 1 << " ";
        }
        cout << endl << endl;
        //    nangles= std::intptr_t (lammps_extract_global(handle, "nangles"));
//    angles = new double[4*nangles];
//    lammps_gather_angles(handle, angles);


//    for (int i=0;i<nangles;++i){
//       for (int j=1;j<4;++j){
//            cout << angles[4*i+j] << " ";
//        }
//        cout << endl;
//    }
        int val;
        cout << endl << endl;
        cout << "Network A Bonds " << a << endl;
        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[a - 1].netCnxs[i] + 1 << " ";
        }
        cout << endl;
        for (int i = 0; i < 3; i++) {
            val = networkA.nodes[a - 1].netCnxs[i] + 1;
            if (val == a) cout << " a  ";
            else if (val == b) cout << " b  ";
            else if (val == c) cout << " c  ";
            else if (val == d) cout << " d  ";
            else if (val == e) cout << " e  ";
            else if (val == f) cout << " f  ";
            else cout << "unk" << endl;
//        cout << networkA.nodes[a-1].netCnxs[i]+1<< " ";
        }
        cout << endl << endl;

        cout << "Network B Bonds" << b << endl;
        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[b - 1].netCnxs[i] + 1 << " ";
        }
        cout << endl;
        for (int i = 0; i < 3; i++) {
            val = networkA.nodes[b - 1].netCnxs[i] + 1;
            if (val == a) cout << " a  ";
            else if (val == b) cout << " b  ";
            else if (val == c) cout << " c  ";
            else if (val == d) cout << " d  ";
            else if (val == e) cout << " e  ";
            else if (val == f) cout << " f  ";
            else cout << "unk" << endl;
//        cout << networkA.nodes[a-1].netCnxs[i]+1<< " ";
        }
        cout << endl << endl;

        cout << "new Lammps Bonds A" << endl;
        cout << b << " " << c << " " << d << endl;
        cout << "b  c  d" << endl;
        cout << endl;
        cout << "new Lammps Bonds B" << endl;
        cout << a << " " << e << " " << f << endl;
        cout << "a  e  f" << endl;
        cout << endl;


        cout << "Network A Angles" << endl;
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                int val0 = networkA.nodes[a - 1].netCnxs[i] + 1;
                int val1 = networkA.nodes[a - 1].netCnxs[j] + 1;
                if (val0 == a) cout << " a  ";
                else if (val0 == b) cout << " b  ";
                else if (val0 == c) cout << " c  ";
                else if (val0 == d) cout << " d  ";
                else if (val0 == e) cout << " e  ";
                else if (val0 == f) cout << " f  ";
                cout << "  a  ";
                if (val1 == a) cout << " a  ";
                else if (val1 == b) cout << " b  ";
                else if (val1 == c) cout << " c  ";
                else if (val1 == d) cout << " d  ";
                else if (val1 == e) cout << " e  ";
                else if (val1 == f) cout << " f  ";
                cout << endl;
//            cout << networkA.nodes[a-1].netCnxs[i] + 1 << " " << a << " " << networkA.nodes[a-1].netCnxs[j] + 1 << endl;
            }
        }
        cout << endl;

        cout << "Lammps Angles A" << endl;
        cout << "c  a  b" << endl;
        cout << "b  a  d" << endl;
        cout << "c  a  d" << endl;
//    cout << c << " " << a << " " << b << endl;
//    cout << d << " " << a << " " << b << endl;
//    cout << d << " " << a << " " << c << endl;
        cout << endl;

        cout << "Network B Angles" << endl;
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                int val0 = networkA.nodes[b - 1].netCnxs[i] + 1;
                int val1 = networkA.nodes[b - 1].netCnxs[j] + 1;
                if (val0 == a) cout << " a  ";
                else if (val0 == b) cout << " b  ";
                else if (val0 == c) cout << " c  ";
                else if (val0 == d) cout << " d  ";
                else if (val0 == e) cout << " e  ";
                else if (val0 == f) cout << " f  ";
                cout << "  b  ";
                if (val1 == a) cout << " a  ";
                else if (val1 == b) cout << " b  ";
                else if (val1 == c) cout << " c  ";
                else if (val1 == d) cout << " d  ";
                else if (val1 == e) cout << " e  ";
                else if (val1 == f) cout << " f  ";
                cout << endl;
//            cout << networkA.nodes[a-1].netCnxs[i] + 1 << " " << a << " " << networkA.nodes[a-1].netCnxs[j] + 1 << endl;
            }
        }
        cout << endl;

        cout << "Lammps Angles B" << endl;
        cout << "a  b  f " << endl;
        cout << "e  b  f " << endl;
        cout << "a  b  e" << endl;
//    cout << a << " " << b << " " << f << endl;
//    cout << e << " " << b << " " << f << endl;
//    cout << a << " " << b << " " << e << endl;
        cout << endl;
    }

    breakBond(a, e, 1);
    breakBond(b, d, 1);
    formBond(a,d, 1);
    formBond(b,e, 1);

    int e1,e11, d1, d11;
    if (networkA.nodes[e-1].netCnxs[0]==b-1) {
        e1=networkA.nodes[e-1].netCnxs[1]+1;
        e11=networkA.nodes[e-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[e-1].netCnxs[1]==b-1) {
        e1=networkA.nodes[e-1].netCnxs[0]+1;
        e11=networkA.nodes[e-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[e-1].netCnxs[2]==b-1) {
        e1=networkA.nodes[e-1].netCnxs[0]+1;
        e11=networkA.nodes[e-1].netCnxs[1]+1;
    }
    else{
        cout << "e not connected to b" << endl;
        sleep(10);
    }
    if (networkA.nodes[d-1].netCnxs[0]==a-1) {
        d1=networkA.nodes[d-1].netCnxs[1]+1;
        d11=networkA.nodes[d-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[d-1].netCnxs[1]==a-1) {
        d1=networkA.nodes[d-1].netCnxs[0]+1;
        d11=networkA.nodes[d-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[d-1].netCnxs[2]==a-1) {
        d1=networkA.nodes[d-1].netCnxs[0]+1;
        d11=networkA.nodes[d-1].netCnxs[1]+1;
    }
    else{
        cout << "d not connected to a" << endl;
        sleep(10);
    }

    breakAngle(c,a,e);
    breakAngle(b,a,e);
    breakAngle(d,b,f);
    breakAngle(d,b,a);

    breakAngle(e1, e, a);
    breakAngle(e11, e, a);
    breakAngle(d1, d, b);
    breakAngle(d11, d, b);

    formAngle(c,a,d);
    formAngle(b,a,d);
    formAngle(f,b,e);
    formAngle(a,b,e);

    formAngle(e1,e,b);
    formAngle(e11,e,b);
    formAngle(d1,d,a);
    formAngle(d11,d,a);


    /*
    string command;

    nbonds0= std::intptr_t (lammps_extract_global(handle, "nbonds"));
    command = "group switch id " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    nbonds1= std::intptr_t (lammps_extract_global(handle, "nbonds"));
    if (nbonds1+1!=nbonds0){
        cout << nbonds1+1 << " " << nbonds0 << endl;
        cout << "ERROR -- bond not found" << endl;
    }

    nbonds0= std::intptr_t (lammps_extract_global(handle, "nbonds"));
    command = "group switch id " + to_string(b) + " " + to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    nbonds1= std::intptr_t (lammps_extract_global(handle, "nbonds"));
    if (nbonds1+1!=nbonds0){
        cout << nbonds1+1 << " " << nbonds0 << endl;
        cout << "ERROR -- bond not found" << endl;
    }

    command = "create_bonds single/bond 1 " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    */

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    nangles0= std::intptr_t (lammps_extract_global(handle, "nangles"));
    command = "group switch id " + to_string(a) + " " + to_string(b) + " " + to_string(c) + " " + to_string(d) + " " + to_string(e) + " " + to_string(f);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");
    nangles1= std::intptr_t (lammps_extract_global(handle, "nangles"));
    cout << "deleted all angles, " << nangles0-nangles1 << " vs expected 6" << endl;
    */

    /*
    cout << "                                                              delete " << c <<" "<< a << " "<< e << endl;
    nangles0= std::intptr_t (lammps_extract_global(handle, "nangles"));
    command = "group switch id " + to_string(c) + " " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");
    nangles1= std::intptr_t (lammps_extract_global(handle, "nangles"));
    if (nangles1+1!=nangles0){
        cout << "ERROR -- angle not found" << endl;
    }


    cout << "                                                              delete " << b <<" "<< a << " "<< e << endl;
    nangles0= std::intptr_t (lammps_extract_global(handle, "nangles"));
    command = "group switch id " + to_string(b) + " " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");
    nangles1= std::intptr_t (lammps_extract_global(handle, "nangles"));
    if (nangles1+1!=nangles0){
        cout << "ERROR -- angle not found" << endl;
    }

    cout << "                                                              delete " << a <<" "<< b << " "<< d << endl;
    command = "group switch id " + to_string(a) + " " + to_string(b) + " " + to_string(d);
    nangles0= std::intptr_t (lammps_extract_global(handle, "nangles"));
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");
    nangles1= std::intptr_t (lammps_extract_global(handle, "nangles"));
    if (nangles1+1!=nangles0){
        cout << "ERROR -- angle not found" << endl;
    }

    cout << "                                                              delete " << d <<" "<< b << " "<< f << endl;
    command = "group switch id " + to_string(d) + " " + to_string(b) + " " + to_string(f);
    nangles0= std::intptr_t (lammps_extract_global(handle, "nangles"));
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");
    nangles1= std::intptr_t (lammps_extract_global(handle, "nangles"));
    if (nangles1+1!=nangles0){
        cout << "ERROR -- angle not found" << endl;
    }

    cout << "                                                              create " << c <<" "<< a << " "<< d << endl;
    command = "create_bonds single/angle 1 " + to_string(c) + " " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());

    cout << "                                                              create " << b <<" "<< a << " "<< d << endl;
    command = "create_bonds single/angle 1 " + to_string(b) + " " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());

    cout << "                                                              create " << a <<" "<< b << " "<< e << endl;
    command = "create_bonds single/angle 1 " + to_string(a) + " " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());

    cout << "                                                              create " << e <<" "<< b << " "<< f << endl;
    command = "create_bonds single/angle 1 " + to_string(f) + " " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    */
}

void LammpsObject::revertGraphene(VecF<int> switchIdsA, Network networkA) {
    if (globalVerbose) cout << "----------------- Revert Graphene" << endl;
    bool verbose=false;

    //unpck parameters
    int a, b, c, d, e, f;
    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];

    a += 1;
    b += 1;
    c += 1;
    d += 1;
    e += 1;
    f += 1;

    /*
    if (a==266 || b==266 || c==266 || d==266 || e==266 || f==266 || a==227 || b==227 || c==227 || d==227 || e==227 || f==227) {
        cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;

        sleep(10);
        cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n INCORPORATED 267 \n\n\n" << endl;
    }
*/
    if (verbose) {
        cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;

        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[a - 1].netCnxs[i] + 1 << " ";
        }
        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[b - 1].netCnxs[i] + 1 << " ";
        }
        cout << endl << endl;
        //    nangles= std::intptr_t (lammps_extract_global(handle, "nangles"));
//    angles = new double[4*nangles];
//    lammps_gather_angles(handle, angles);


//    for (int i=0;i<nangles;++i){
//       for (int j=1;j<4;++j){
//            cout << angles[4*i+j] << " ";
//        }
//        cout << endl;
//    }



        int val;
        cout << endl << endl;
        cout << "Network A Bonds " << a << endl;
        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[a - 1].netCnxs[i] + 1 << " ";
        }
        cout << endl;
        for (int i = 0; i < 3; i++) {
            val = networkA.nodes[a - 1].netCnxs[i] + 1;
            if (val == a) cout << " a  ";
            else if (val == b) cout << " b  ";
            else if (val == c) cout << " c  ";
            else if (val == d) cout << " d  ";
            else if (val == e) cout << " e  ";
            else if (val == f) cout << " f  ";
            else cout << "unk" << endl;
//        cout << networkA.nodes[a-1].netCnxs[i]+1<< " ";
        }
        cout << endl << endl;

        cout << "Network B Bonds" << b << endl;
        for (int i = 0; i < 3; i++) {
            cout << networkA.nodes[b - 1].netCnxs[i] + 1 << " ";
        }
        cout << endl;
        for (int i = 0; i < 3; i++) {
            val = networkA.nodes[b - 1].netCnxs[i] + 1;
            if (val == a) cout << " a  ";
            else if (val == b) cout << " b  ";
            else if (val == c) cout << " c  ";
            else if (val == d) cout << " d  ";
            else if (val == e) cout << " e  ";
            else if (val == f) cout << " f  ";
            else cout << "unk" << endl;
//        cout << networkA.nodes[a-1].netCnxs[i]+1<< " ";
        }
        cout << endl << endl;

        cout << "new Lammps Bonds A" << endl;
        cout << b << " " << c << " " << e << endl;
        cout << "b  c  e" << endl;
        cout << endl;
        cout << "new Lammps Bonds B" << endl;
        cout << a << " " << d << " " << f << endl;
        cout << "a  d  f" << endl;
        cout << endl;


        cout << "Network A Angles" << endl;
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                int val0 = networkA.nodes[a - 1].netCnxs[i] + 1;
                int val1 = networkA.nodes[a - 1].netCnxs[j] + 1;
                if (val0 == a) cout << " a  ";
                else if (val0 == b) cout << " b  ";
                else if (val0 == c) cout << " c  ";
                else if (val0 == d) cout << " d  ";
                else if (val0 == e) cout << " e  ";
                else if (val0 == f) cout << " f  ";
                cout << "  a  ";
                if (val1 == a) cout << " a  ";
                else if (val1 == b) cout << " b  ";
                else if (val1 == c) cout << " c  ";
                else if (val1 == d) cout << " d  ";
                else if (val1 == e) cout << " e  ";
                else if (val1 == f) cout << " f  ";
                cout << endl;
//            cout << networkA.nodes[a-1].netCnxs[i] + 1 << " " << a << " " << networkA.nodes[a-1].netCnxs[j] + 1 << endl;
            }
        }
        cout << endl;

        cout << "Lammps Angles A" << endl;
        cout << "c  a  e" << endl;
        cout << "b  a  e" << endl;
        cout << "c  a  b" << endl;
//    cout << c << " " << a << " " << b << endl;
//    cout << d << " " << a << " " << b << endl;
//    cout << d << " " << a << " " << c << endl;
        cout << endl;

        cout << "Network B Angles" << endl;
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                int val0 = networkA.nodes[b - 1].netCnxs[i] + 1;
                int val1 = networkA.nodes[b - 1].netCnxs[j] + 1;
                if (val0 == a) cout << " a  ";
                else if (val0 == b) cout << " b  ";
                else if (val0 == c) cout << " c  ";
                else if (val0 == d) cout << " d  ";
                else if (val0 == e) cout << " e  ";
                else if (val0 == f) cout << " f  ";
                cout << "  b  ";
                if (val1 == a) cout << " a  ";
                else if (val1 == b) cout << " b  ";
                else if (val1 == c) cout << " c  ";
                else if (val1 == d) cout << " d  ";
                else if (val1 == e) cout << " e  ";
                else if (val1 == f) cout << " f  ";
                cout << endl;
//            cout << networkA.nodes[a-1].netCnxs[i] + 1 << " " << a << " " << networkA.nodes[a-1].netCnxs[j] + 1 << endl;
            }
        }
        cout << endl;

        cout << "Lammps Angles B" << endl;
        cout << "a  b  f " << endl;
        cout << "d  b  f " << endl;
        cout << "a  b  d" << endl;
//    cout << a << " " << b << " " << f << endl;
//    cout << e << " " << b << " " << f << endl;
//    cout << a << " " << b << " " << e << endl;
        cout << endl;
    }

    breakBond(a,d,1);
    breakBond(b,e,1);
    formBond(a,e,1);
    formBond(b,d,1);

    int e1,e11, d1, d11;
    if (networkA.nodes[e-1].netCnxs[0]==a-1) {
        e1=networkA.nodes[e-1].netCnxs[1]+1;
        e11=networkA.nodes[e-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[e-1].netCnxs[1]==a-1) {
        e1=networkA.nodes[e-1].netCnxs[0]+1;
        e11=networkA.nodes[e-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[e-1].netCnxs[2]==a-1) {
        e1=networkA.nodes[e-1].netCnxs[0]+1;
        e11=networkA.nodes[e-1].netCnxs[1]+1;
    }
    else{
        cout << "e not connected to a" << endl;
        sleep(10);
    }
    if (networkA.nodes[d-1].netCnxs[0]==b-1) {
        d1=networkA.nodes[d-1].netCnxs[1]+1;
        d11=networkA.nodes[d-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[d-1].netCnxs[1]==b-1) {
        d1=networkA.nodes[d-1].netCnxs[0]+1;
        d11=networkA.nodes[d-1].netCnxs[2]+1;
    }
    else if (networkA.nodes[d-1].netCnxs[2]==b-1) {
        d1=networkA.nodes[d-1].netCnxs[0]+1;
        d11=networkA.nodes[d-1].netCnxs[1]+1;
    }
    else{
        cout << "d not connected to b" << endl;
        sleep(10);
    }



    breakAngle(c,a,d);
    breakAngle(b,a,d);
    breakAngle(f,b,e);
    breakAngle(a,b,e);

    breakAngle(e1,e,b);
    breakAngle(e11,e,b);
    breakAngle(d1,d,a);
    breakAngle(d11,d,a);

    formAngle(c,a,e);
    formAngle(b,a,e);
    formAngle(d,b,f);
    formAngle(d,b,a);

    formAngle(e1, e, a);
    formAngle(e11, e, a);
    formAngle(d1, d, b);
    formAngle(d11, d, b);

    /*
    string command;
    command = "group switch id " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b) + " " + to_string(d);
    lammps_command(handle, command.c_str());


    command = "group switch id " + to_string(c) + " " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(a) + " " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(f) + " " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/angle 1 " + to_string(c) + " " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/angle 1 " + to_string(b) + " " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/angle 1 " + to_string(a) + " " + to_string(b) + " " + to_string(d);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/angle 1 " + to_string(d) + " " + to_string(b) + " " + to_string(f);
    lammps_command(handle, command.c_str());
     */
}


void LammpsObject::switchTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT, Network networkT){
    if (globalVerbose) cout << "----------------- Switch Triangle Raft" << endl;

    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    int alpha, beta, gamma, delta, eta;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    alpha =     switchIdsT[6];
    beta =      switchIdsT[7];
    gamma =     switchIdsT[8];
    delta =     switchIdsT[9];
    eta =       switchIdsT[10];



    a +=1;
    b +=1;

    alpha +=1;
    beta  +=1;
    gamma +=1;
    delta +=1;
    eta   +=1;


    //Corrections to
    // -- Actually remove the bonds
    // -- Correct bond types
    
    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */


    //cout << " a - delta" << endl;
    //breakBond(a, alpha);

    breakBond(a,gamma,2);
    breakBond(b, delta,2);
    breakBond(gamma, beta,1);
    breakBond(delta, eta,1);

    /*
    string command;
    command = "group switch id " + to_string(a) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(delta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma) + " " + to_string(beta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta) + " " + to_string(eta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */

    formBond(a,delta,2);
    formBond(b, gamma,2);
    formBond(beta, delta,1);
    formBond(eta, gamma,1);
    /*
    command = "create_bonds single/bond 2 " + to_string(a) + " " + to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 2 " + to_string(b) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta) + " " + to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());
    */

}

void LammpsObject::revertTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT){
    if (globalVerbose) cout << "----------------- Revert Triangle Raft" << endl;
    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    int alpha, beta, gamma, delta, eta;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    alpha =     switchIdsT[6];
    beta =      switchIdsT[7];
    gamma =     switchIdsT[8];
    delta =     switchIdsT[9];
    eta =       switchIdsT[10];


    a +=1;
    b +=1;

    alpha +=1;
    beta  +=1;
    gamma +=1;
    delta +=1;
    eta   +=1;


    breakBond(a,delta,2);
    breakBond(b, gamma,2);
    breakBond(beta, delta,1);
    breakBond(eta, gamma,1);
    /*
    *
    *      A x gamma       A - delta
    *      B x delta       B - gamma
    *
    *      gamma x beta    gamma - eta
    *      beta  x gammma  beta - delta
    *      delta x eta     delta - beta
    *      eta   x delta   eta - gamma
    *
    */


    /*
    string command;
    command = "group switch id " + to_string(a) + " " + to_string(delta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma) + " " + to_string(eta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta) + " " + to_string(beta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */

    formBond(a,gamma,2);
    formBond(b, delta,2);
    formBond(gamma, beta,1);
    formBond(delta, eta,1);
    /*
    command = "create_bonds single/bond 2 " + to_string(a) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 2 " + to_string(b) + " " + to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta) + " " + to_string(delta);
    lammps_command(handle, command.c_str());
     */
}


void LammpsObject::switchBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT){

    if (globalVerbose) cout << "----------------- Switch Bilayer" << endl;
    bool verbose=false;
    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    int alpha, beta, gamma, delta, eta;

    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    alpha =     switchIdsT[6];
    beta =      switchIdsT[7];
    gamma =     switchIdsT[8];
    delta =     switchIdsT[9];
    eta =       switchIdsT[10];


    /*
  *
  *      A x gamma       A - delta
  *      B x delta       B - gamma
  *
  *      gamma x beta    gamma - eta
  *      beta  x gammma  beta - delta
  *      delta x eta     delta - beta
  *      eta   x delta   eta - gamma
  *
  */

    int nSi, nSi2, nO;
//    natoms = std::intptr_t (lammps_get_natoms(handle));

//    cout << "natoms standard def " << natoms << " " << natoms/19164 << endl;
//    double natoms2 = lammps_get_natoms(handle);

    
//    double *natoms3 = intptr_t (lammps_extract_compute(handle, "1 all count/type atom", 0, 1));
//    double *natoms5 = (double *)lammps_extract_compute(handle, "1 all count/type atom", 0, 1);
//    double *natoms3 = (double *)lammps_extract_compute(handle, "atomcount", 0, 1);

//    cout << "pointer to atomcount[0] " << natoms3[0] << endl;
//    int natoms4 = sizeof(natoms3)/sizeof(natoms3[0]);
//    cout << "size of array calculation : " << natoms4 << endl;



//    double *natoms4 = (double *)lammps_extract_compute(handle, "1 all count/type {atom}", 0, 0);

//    cout << natoms3 << endl;
//   (int)(lammps_get_natoms(handle)+0.5)
    nSi = (int)(round(natoms/3)+0.5);
    nSi2= (int)(round(nSi/2)+0.5);
    nO =  natoms - nSi;
    

//    cout << natoms << " " << nSi << " " << nO << endl;
//    cout << natoms2 << " " << natoms2/3 << " " << round(natoms2/3) << endl;
//    cout << nO << " " << 2*natoms2/3 << " " << round(2*natoms2/3) << endl;
//    sleep(100);


    int a_prime, a_prime_prime, b_prime, b_prime_prime;
    int a_prime_prime_prime, b_prime_prime_prime;

    int beta_prime, beta_prime_prime, gamma_prime, gamma_prime_prime;
    int eta_prime, eta_prime_prime, delta_prime, delta_prime_prime;

    a_prime = 2*a+1;
    a_prime_prime = 2*a+2;
    a_prime_prime_prime = nSi+a+1;

    b_prime = 2*b+1;
    b_prime_prime = 2*b+2;
    b_prime_prime_prime = nSi+b+1;

    beta_prime = 3*nSi2 + (beta - nSi2)*2+1;
    beta_prime_prime = beta_prime +1;
    gamma_prime = 3*nSi2 + (gamma - nSi2)*2+1;
    gamma_prime_prime = gamma_prime +1;
    eta_prime = 3*nSi2 + (eta - nSi2)*2+1;
    eta_prime_prime = eta_prime +1;
    delta_prime = 3*nSi2 + (delta - nSi2)*2+1;
    delta_prime_prime = delta_prime +1;

    if (verbose){
        cout << a << " " << b << endl;
        cout << a_prime << " " << a_prime_prime << " " << a_prime_prime_prime << endl;
        cout << b_prime << " " << b_prime_prime << " " << b_prime_prime_prime << endl;

        cout << endl;
        cout << beta << " " << gamma << " " << eta << " " << delta << endl;
        cout << beta_prime << " " << beta_prime_prime << endl;
        cout << gamma_prime << " " << gamma_prime_prime << endl;
        cout << eta_prime << " " << eta_prime_prime << endl;
        cout << delta_prime << " " <<delta_prime_prime << endl;



    }


    //      Top layer
    //cout << "Top Layer" << endl;

    breakBond(a_prime, gamma_prime,2);
    breakBond(b_prime, delta_prime,2);
    breakBond(beta_prime, gamma_prime,1);
    breakBond(eta_prime, delta_prime,1);
    /*

    string command;
    command = "group switch id " + to_string(a_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma_prime) + " " + to_string(beta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta_prime) + " " + to_string(eta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */
    formBond(a_prime, delta_prime,2);
    formBond(b_prime, gamma_prime,2);
    formBond(delta_prime, beta_prime,1);
    formBond(gamma_prime, eta_prime,1);
    /*
    command = "create_bonds single/bond 2 " + to_string(a_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 2 " + to_string(b_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    */

    // Bottom layer
    //cout << "Bottom Layer" << endl;

    breakBond(a_prime_prime, gamma_prime_prime,2);
    breakBond(b_prime_prime, delta_prime_prime,2);
    breakBond(beta_prime_prime, gamma_prime_prime,1);
    breakBond(eta_prime_prime, delta_prime_prime,1);
    /*
    string command;
    command = "group switch id " + to_string(a_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma_prime_prime) + " " + to_string(beta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta_prime_prime) + " " + to_string(eta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */
    formBond(a_prime_prime, delta_prime_prime,2);
    formBond(b_prime_prime, gamma_prime_prime,2);
    formBond(beta_prime_prime, delta_prime_prime,1);
    formBond(eta_prime_prime, gamma_prime_prime, 1);
    /*
    command = "create_bonds single/bond 2 " + to_string(a_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 2 " + to_string(b_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    */
    // Bridge to Bottom Layer
    //cout << "Bridge to Bottom Layer" << endl;
    breakBond(a_prime_prime_prime, gamma_prime_prime,1);
    breakBond(b_prime_prime_prime, delta_prime_prime,1);
    /*
    command = "group switch id " + to_string(a_prime_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */

    formBond(a_prime_prime_prime, delta_prime_prime,1);
    formBond(b_prime_prime_prime, gamma_prime_prime,1);
    /*
    command = "create_bonds single/bond 1 " + to_string(a_prime_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    */

    // Bridge to Top Layer
    //cout << "Bridge to Top Layer" << endl;

    breakBond(a_prime_prime_prime, gamma_prime,1);
    breakBond(b_prime_prime_prime, delta_prime,1);
    /*
    command = "group switch id " + to_string(a_prime_prime_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */
    formBond(a_prime_prime_prime, delta_prime,1);
    formBond(b_prime_prime_prime, gamma_prime,1 );
    /*
    command = "create_bonds single/bond 1 " + to_string(a_prime_prime_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime_prime_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    */



}

void LammpsObject::revertBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT){

    if (globalVerbose) cout << "----------------- Revert Bilayer" << endl;
    bool verbose=false;
    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    int alpha, beta, gamma, delta, eta;

    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    alpha =     switchIdsT[6];
    beta =      switchIdsT[7];
    gamma =     switchIdsT[8];
    delta =     switchIdsT[9];
    eta =       switchIdsT[10];

    /*
  *
  *      A x gamma       A - delta
  *      B x delta       B - gamma
  *
  *      gamma x beta    gamma - eta
  *      beta  x gammma  beta - delta
  *      delta x eta     delta - beta
  *      eta   x delta   eta - gamma
  *
  */

    int nSi, nSi2, nO; 

    nSi = (int)(round(natoms/3)+0.5);
    nSi2= (int)(round(nSi/2)+0.5);
    nO =  natoms - nSi;

    int a_prime, a_prime_prime, b_prime, b_prime_prime;
    int a_prime_prime_prime, b_prime_prime_prime;

    int beta_prime, beta_prime_prime, gamma_prime, gamma_prime_prime;
    int eta_prime, eta_prime_prime, delta_prime, delta_prime_prime;

    a_prime = 2*a+1;
    a_prime_prime = 2*a+2;
    a_prime_prime_prime = nSi+a+1;

    b_prime = 2*b+1;
    b_prime_prime = 2*b+2;
    b_prime_prime_prime = nSi+b+1;

    beta_prime = 3*nSi2 + (beta - nSi2)*2+1;
    beta_prime_prime = beta_prime +1;
    gamma_prime = 3*nSi2 + (gamma - nSi2)*2+1;
    gamma_prime_prime = gamma_prime +1;
    eta_prime = 3*nSi2 + (eta - nSi2)*2+1;
    eta_prime_prime = eta_prime +1;
    delta_prime = 3*nSi2 + (delta - nSi2)*2+1;
    delta_prime_prime = delta_prime +1;

    /*
     * a_prime = nO+2*a+1;
    a_prime_prime = nO+2*a+2;
    a_prime_prime_prime = a+1;

    b_prime = nO+2*b+1;
    b_prime_prime = nO+2*b+2;
    b_prime_prime_prime = b+1;

    beta_prime = 3*nSi2 + (beta - nSi2)*2+1;
    beta_prime_prime = beta_prime +1;
    gamma_prime = 3*nSi2 + (gamma - nSi2)*2+1;
    gamma_prime_prime = gamma_prime +1;
    eta_prime = 3*nSi2 + (eta - nSi2)*2+1;
    eta_prime_prime = eta_prime +1;
    delta_prime = 3*nSi2 + (delta - nSi2)*2+1;
    delta_prime_prime = delta_prime +1;
    */

    if (verbose){
        cout << a << " " << b << endl;
        cout << a_prime << " " << a_prime_prime << " " << a_prime_prime_prime << endl;
        cout << b_prime << " " << b_prime_prime << " " << b_prime_prime_prime << endl;

        cout << endl;
        cout << beta << " " << gamma << " " << eta << " " << delta << endl;
        cout << beta_prime << " " << beta_prime_prime << endl;
        cout << gamma_prime << " " << gamma_prime_prime << endl;
        cout << eta_prime << " " << eta_prime_prime << endl;
        cout << delta_prime << " " <<delta_prime_prime << endl;



    }



    /*
    a_prime = 2*a;
    a_prime_prime = 2*a+1;
    a_prime_prime_prime = 3*a;

    b_prime = 2*b;
    b_prime_prime = 2*b+1;
    b_prime_prime_prime = 3*b;

    beta_prime = 3*nSi + (beta - nSi)*2;
    beta_prime_prime = beta_prime +1;
    gamma_prime = 3*nSi + (gamma - nSi)*2;
    gamma_prime_prime = gamma_prime +1;
    eta_prime = 3*nSi + (eta - nSi)*2;
    eta_prime_prime = eta_prime +1;
    delta_prime = 3*nSi + (delta - nSi)*2;
    delta_prime_prime = delta_prime +1;
    */

    //      Top layer
    //cout << "Top Layer" << endl;
    string command;

    breakBond(a_prime, delta_prime,2);
    breakBond(b_prime, gamma_prime,2);
    breakBond(delta_prime, beta_prime,1);
    breakBond(gamma_prime, eta_prime,1);
    /*
    command = "group switch id " + to_string(a_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta_prime) + " " + to_string(beta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma_prime) + " " + to_string(eta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */
    formBond(a_prime, gamma_prime,2);
    formBond(b_prime, delta_prime,2);
    formBond(beta_prime, gamma_prime,1);
    formBond(eta_prime, delta_prime,1);
    /*
    command = "create_bonds single/bond 2 " + to_string(a_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 2 " + to_string(b_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    */



    // Bottom layer
    //cout << " Bottom Layer" << endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    breakBond(a_prime_prime, delta_prime_prime,2);
    breakBond(b_prime_prime, gamma_prime_prime,2);
    breakBond(delta_prime_prime, beta_prime_prime,1);
    breakBond(gamma_prime_prime, eta_prime_prime,1);
    /*
    command = "group switch id " + to_string(a_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 2 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta_prime_prime) + " " + to_string(beta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma_prime_prime) + " " + to_string(eta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */
    formBond(a_prime_prime, gamma_prime_prime,2);
    formBond(b_prime_prime, delta_prime_prime,2);
    formBond(beta_prime_prime, gamma_prime_prime,1);
    formBond(eta_prime_prime, delta_prime_prime, 1);
    /*
    command = "create_bonds single/bond 2 " + to_string(a_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 2 " + to_string(b_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    */

    // Bridge to Top Layer

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //cout << "Bridge to Top Layer" << endl;
    breakBond(a_prime_prime_prime, delta_prime,1);
    breakBond(b_prime_prime_prime, gamma_prime,1);

    /*
    command = "group switch id " + to_string(a_prime_prime_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */

    formBond(a_prime_prime_prime, gamma_prime,1);
    formBond(b_prime_prime_prime, delta_prime, 1);
    /*
    command = "create_bonds single/bond 1 " + to_string(a_prime_prime_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime_prime_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    */




    // Bridge to Bottom Layer
    //cout << "Bridge to Bottom Layer" << endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    breakBond(a_prime_prime_prime, delta_prime_prime,1);
    breakBond(b_prime_prime_prime, gamma_prime_prime,1);
    /*
    command = "group switch id " + to_string(a_prime_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch bond 1 remove");
    lammps_command(handle, "group switch delete");
    */
    formBond(a_prime_prime_prime, gamma_prime_prime,1);
    formBond(b_prime_prime_prime, delta_prime_prime, 1);
    /*
    command = "create_bonds single/bond 1 " + to_string(a_prime_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    */


}



void LammpsObject::switchBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT){
    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    int alpha, beta, gamma, delta, eta;

    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    alpha =     switchIdsT[6];
    beta =      switchIdsT[7];
    gamma =     switchIdsT[8];
    delta =     switchIdsT[9];
    eta =       switchIdsT[10];

   /* Switch connectivities in lattice and dual
    * 3-3 coordination connection
    * a,b,c,d,e,f are nodes in lattice A
    * u,v,w,x are nodes in lattice B
    *  E      F            V
    *   \    /           / | \
    *   A---B           W  |  X
    *  /     \           \ | /
    * C       D            U
    *
    *       A x E       A - D
    *       B x D       B - E
    *       D x B       D - A
    *       E x A       E - B
    *
    *       Break A - E and B - D
    *       Form  A - D and B - E
    *
    */

    string command;
    command = "group switch id " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());


    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */
    command = "group switch id " + to_string(a) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(delta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma) + " " + to_string(beta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta) + " " + to_string(eta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");


    command = "create_bonds single/bond 1 " + to_string(a) + " " + to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta) + " " + to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());


   /*
    *
    *      A x gamma       A - delta
    *      B x delta       B - gamma
    *
    *      gamma x beta    gamma - eta
    *      beta  x gammma  beta - delta
    *      delta x eta     delta - beta
    *      eta   x delta   eta - gamma
    *
    */

    int nSi, nO;
    nSi = round(natoms/3);
    nO =  round(2*natoms/3);

    int a_prime, a_prime_prime, b_prime, b_prime_prime;
    int a_prime_prime_prime, b_prime_prime_prime;

    int beta_prime, beta_prime_prime, gamma_prime, gamma_prime_prime;
    int eta_prime, eta_prime_prime, delta_prime, delta_prime_prime;

    a_prime = 2*a;
    a_prime_prime = 2*a+1;
    a_prime_prime_prime = 3*a;
    
    b_prime = 2*b;
    b_prime_prime = 2*b+1;
    b_prime_prime_prime = 3*b;
    
    beta_prime = 3*nSi + (beta - nSi)*2;
    beta_prime_prime = beta_prime +1;
    gamma_prime = 3*nSi + (gamma - nSi)*2;
    gamma_prime_prime = gamma_prime +1;
    eta_prime = 3*nSi + (eta - nSi)*2;
    eta_prime_prime = eta_prime +1;
    delta_prime = 3*nSi + (delta - nSi)*2;
    delta_prime_prime = delta_prime +1;

    //      Top layer

    command = "group switch id " + to_string(a_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma_prime) + " " + to_string(beta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta_prime) + " " + to_string(eta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");


    command = "create_bonds single/bond 1 " + to_string(a_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    // Bottom layer
    
    command = "group switch id " + to_string(a_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma_prime_prime) + " " + to_string(beta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta_prime_prime) + " " + to_string(eta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");


    command = "create_bonds single/bond 1 " + to_string(a_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    // Bridge to Top Layer

    command = "group switch id " + to_string(a_prime_prime_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + to_string(a_prime_prime_prime) + " " + to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime_prime_prime) + " " + to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    // Bridge to Bottom Layer

    command = "group switch id " + to_string(a_prime_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b_prime_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + to_string(a_prime_prime_prime) + " " + to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b_prime_prime_prime) + " " + to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());





}

void LammpsObject::revertBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT){
    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    int alpha, beta, gamma, delta, eta;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    alpha =     switchIdsT[6];
    beta =      switchIdsT[7];
    gamma =     switchIdsT[8];
    delta =     switchIdsT[9];
    eta =       switchIdsT[10];
    /* Switch connectivities in lattice and dual
     * 3-3 coordination connection
     * a,b,c,d,e,f are nodes in lattice A
     * u,v,w,x are nodes in lattice B
     *  E      F            V
     *   \    /           / | \
     *   A---B           W  |  X
     *  /     \           \ | /
     * C       D            U
     *
     *       A x E       A - D
     *       B x D       B - E
     *       D x B       D - A
     *       E x A       E - B
     *
     *       Break A - E and B - D
     *       Form  A - D and B - E
     *
     */

    string command;
    command = "group switch id " + to_string(a) + " " + to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + to_string(a) + " " + to_string(e);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b) + " " + to_string(d);
    lammps_command(handle, command.c_str());


    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */
    command = "group switch id " + to_string(a) + " " + to_string(delta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(b) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(gamma) + " " + to_string(eta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + to_string(delta) + " " + to_string(beta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");


    command = "create_bonds single/bond 1 " + to_string(a) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(b) + " " + to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(beta) + " " + to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + to_string(eta) + " " + to_string(delta);
    lammps_command(handle, command.c_str());
}

VecF<int> LammpsObject::GlobalPotentialMinimisation(){
/*
    lammps_command(handle, "min_style               cg\n"
                           "minimize        1.0e-6 0.0 1000000 10000000\n"
                           "\n"
                           "min_style               sd\n"
                           "minimize        1.0e-6 0.0 1000000 10000000\n");
*/

    lammps_command(handle, "minimize 1.0e-6 0.0 1000000 10000000");
//    cout << "minimize" << endl;


    VecF<int> optstatus(2);
    optstatus[0]=0;
    optstatus[1]=10;

    return optstatus;
}

double LammpsObject::GlobalPotentialEnergy() {
//    lammps_command(handle, "min_style               cg\n"
//                           "minimize        1.0e-6 0.0 1000000 10000000\n");
//    double *dptr = (double *) lammps_extract_global(handle, )
//    return lammps_extract_global()

//    lammps_command(handle,"run 100 pre no post no");
//    printf("PE = %g\nKE = %g\n",
//           lammps_get_thermo(handle,"pe"),
//           lammps_get_thermo(handle,"ke"));

    return lammps_get_thermo(handle, "pe");
}

/*
void LammpsObject::groupAtoms(VecR<int> switch_atoms){
    // feed vector of atoms to change

    string command;

    for (int i=0; i<switch_atoms.n;++i){
        command = "group switch id "+ to_string(switch_atoms[i]);
        lammps_command(handle, command.c_str());
    }
    // group atoms to delete bonds
    // group ID style args
    //lammps_command(handle, command);

    // delete all bonds and angles
    lammps_command(handle, "delete switch bond 1 remove");
    lammps_command(handle, "delete switch bond 2 remove");
    lammps_command(handle, "delete switch angle 1 remove");

    lammps_command(handle, "delete switch angle 1 remove");

    // recreate bonds
    // loop over atoms in linked network, if both atoms in 'switch' and id1>id0, create the bond
    for (int id0 =0; id0<network.nodes.n;++j){
        for (int i = 0; i < network.nodes[id0].netCnxs.n; ++i) {
            id1 = network.nodes[id0].netCnxs[i];
            id2 = network.nodes[id0].netCnxs[(i + 1) % cnd];
            //bonds
            if (id0 < id1) {
                if (network.nodes.n>networkA.nodes.n){
                    //prevent double counting
                    //if (id0>networkA.nodes.n){cout<<"large"<<endl;}
                    //if (id1>networkA.nodes.n){cout<<"large"<<endl;}
                    if (id0>=networkA.nodes.n && id1>=networkA.nodes.n){
                        bonds.addValue(id0);
                        bonds.addValue(id1);
                        bondParams.addValue(bk);
                        bondParams.addValue(br0/2);
                        sio_bond_count +=1;
//                    cout << id0 << " " << id1 << " " << bk << " " << br0/2 << endl;
                    }
                    else{
                        bonds.addValue(id0);
                        bonds.addValue(id1);
                        bondParams.addValue(bk);
                        bondParams.addValue(br/2);
                        o_o_bond_count +=1;
//                    cout << id0 << " " << id1 << " " << bk << " " << br/2 << endl;

                    }
                }
                else{
                    bonds.addValue(id0);
                    bonds.addValue(id1);
                    bondParams.addValue(bk);
                    bondParams.addValue(br);
//                cout << id0 << " " << id1 << " " << bk << " " << br << endl;

                }


            }




            for (int i=0; i<;++i){

            }

            // recreate angles

            // ungroup
            lammps_command(handle, "group switch delete");

        }
    }
    }
*/


double LammpsObject::probeE(){
    return lammps_get_thermo(handle,"etotal");
}
/*
int main(int argc, char **argv)
{
//    void *handle;
//    int version;


    / delete LAMMPS instance and shut down MPI /

}
*/

/*
void LammpsObject::write(string prefix, string filename, int dim){

    // This will be similar for all systems EXCEPT periodic boundary conditions
    ifstream inputFile(prefix+"_A_aux.dat", ios::in)

    //auxilary information
    ofstream auxFile(prefix + "_" + filename +"_aux.dat", ios::in | ios::trunc);
    auxFile << fixed << showpoint << setprecision(1);

    string skip;
    getline(inputFile,skip);                        // n_nodes
    auxFile << setw(10) << left << natoms <<endl;

    string line;
    int maxdualCnxs, maxnetCnxs;
    getline(inputFile, line);
    maxdualCnxs = line.erase(0,10);

    if (filename=="BL")                    maxnetCnxs=8;
    else if (filename=="TR")               maxnetCnxs=6;
    else                                   maxnetCnxs=3;

    auxFile << setw(10) << left << maxnetCnxs << setw(10) << left << maxdualCnxs << endl;

    string geometryCode;
    getline(inputFile, geometryCode);
    auxFile << setw(10) << left << geometryCode <<endl;

    VecF<double> pb, rpb;
    pb[0] = pbx();
    pb[1] = pby();
    rpb[0] = 1/pb[0];
    rpb[1] = 1/pb[1];
    if (dim>2) {
        pb[2] = pbz();
        rpb[2] = 1/pb[2];
    }


    getline(inputFile, skip);
    getline(inputFile, skip);
    auxFile << fixed << showpoint << setprecision(6);
    for(int i=0; i<pb.n; ++i) auxFile << setw(20) << left << pb[i];
    auxFile<<endl;
    for(int i=0; i<rpb.n; ++i) auxFile << setw(20) << left << rpb[i];
    auxFile <<endl;
    auxFile.close();

    double* write_crds = fetchCrds(dim);
    ofstream crdFile(prefix + "_" + filename + "_crds.dat", ios::in | ios::trunc);
    crdFile << fixed << showpoint << setprecision(6);
    for (int i=0;i<natoms;++i){
        crdFile << setw(20) << left << write_crds[3*i];
        crdFile << setw(20) << left << write_crds[3*i+1];
        crdFile << setw(20) << left << write_crds[3*i+2];
        crdFile << endl;
    }
    crdFile.close()

    ofstream netFile(prefix + "_net.dat", ios::in | ios::trunc);
    netFile << fixed << showpoint << setprecision(1);
    VecF<int> netCnxs(natoms*maxnetCnxs);


    //network connections
    ofstream netFile(prefix + "_net.dat", ios::in | ios::trunc);
    netFile << fixed << showpoint << setprecision(1);
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].netCnxs.n; ++j){
            netFile<< setw(20) << left << nodes[i].netCnxs[j];
        }
        netFile<<endl;
    }
    netFile.close();

}
*/