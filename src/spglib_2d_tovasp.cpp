#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <regex>
#include "spglib.h"

using namespace std;


// Read all lines in the file to a vector
void readlines(ifstream &rfile, vector<string> &buffer){
  string line;
  while (getline(rfile, line)){
    buffer.push_back(line);
  }
}


// Check if the POSCAR file is opened
int poscar_is_not_open(ifstream &target, string filename){
  if (!target.is_open()){
    cerr << "[error] Cannot open the file `./" << filename << "`!" << endl;
    return 1;
  }
  return 0;
}


// Check if the content os POSCAR is valid
int poscar_is_not_valid(vector<string> &buffer){
  string coord_type_line = buffer[7];
  string match_pattern = "Direct";
  // 'Direct' Char not found 
  if (coord_type_line.find(match_pattern) == string::npos) {
    cerr << "[error] The 8th line of POSCAR must be `Direct`!" << endl;
    cerr << "[error] The 6th line of POSCAR must be element symbols" << endl;
    return 1;
  }
  return 0;
}


// Read in header informations, 
//  which include file info. and scaling factor
void readposcar_header(vector<string> &buffer, string &file_info, 
                       string &scale_factor){
  file_info = buffer[0];
  scale_factor = buffer[1];
}


// Read in lattice info.s
void readposcar_lattice(vector<string> &buffer, double lattice_vec[3][3]){
  // - Lattice vector a1,a2,a3
  istringstream buffer_a1(buffer[2]);
  istringstream buffer_a2(buffer[3]);
  istringstream buffer_a3(buffer[4]);
  double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
  // - Split the string of lattice vector
  buffer_a1 >> a1x >> a1y >> a1z;
  buffer_a2 >> a2x >> a2y >> a2z;
  buffer_a3 >> a3x >> a3y >> a3z;
  // - Write lattice vector to lattice matrix
  lattice_vec[0][0] = a1x;
  lattice_vec[0][1] = a2x;
  lattice_vec[0][2] = a3x;
  lattice_vec[1][0] = a1y;
  lattice_vec[1][1] = a2y;
  lattice_vec[1][2] = a3y;
  lattice_vec[2][0] = a1z;
  lattice_vec[2][1] = a2z;
  lattice_vec[2][2] = a3z;
}


// Calculate the total atom number, use in `readposcar_atom_type()`
int calc_total_atom_num(string buffer_elems){
  istringstream iss_buffer_elems(buffer_elems);
  int curr_elem_atom_num;
  int total_atoms_num = 0;
  while (iss_buffer_elems >> curr_elem_atom_num){
    total_atoms_num += curr_elem_atom_num;
  }
  return total_atoms_num;
}

// Make a clear string, del all unnecessary char
void make_clear_string(string &target, string &reult){
  // Clear all extra space in the string
  reult = regex_replace(target, regex("^ +| +$|( ) +"), "$1");
}

// Decided the atom type list, such as, Fe3O4 is 1 1 1 2 2 2 2.
void decided_atom_type_list(string &buffer_elems, int* type_list){
  int curr_elem_atom_num;
  int i_curr_atom;
  int i_curr_elem = 0;
  int tag_curr_elem = 1;
  istringstream iss_buffer_elems(buffer_elems);
  // -- For each element
  while (iss_buffer_elems >> curr_elem_atom_num){
    // -- For each atom of the same element
    for (int i_same_elem=0; i_same_elem<curr_elem_atom_num; ++i_same_elem){
      i_curr_atom = i_curr_elem + i_same_elem;
      type_list[i_curr_atom] = tag_curr_elem;
    }
    i_curr_elem += curr_elem_atom_num;
    ++tag_curr_elem;
  }
}

// Read in the atom type 
//   E.g., Fe2O3 is 1 1 2 2 2
void readposcar_atom_type(vector<string> &buffer, int* &type_list, 
                          int &atoms_num){
  // - Read the 7th line for atoms number of each element
  string buffer_elems;
  make_clear_string(buffer[6], buffer_elems); // Del. all uesless blank
  // - Calculate the total atom number
  atoms_num = calc_total_atom_num(buffer_elems);
  // - Write the type list
  type_list = new int[atoms_num];
  decided_atom_type_list(buffer_elems, type_list);
}


// Get the x-th atom's line number in the POSCAR
int get_xth_atom_line_num_in_poscar(int i_atom){
  const int ATOM_COORD_LINE = 8; // line 9 - 1
  return ATOM_COORD_LINE + i_atom;
}

// Read string to istringsteam vals
void str_into_istringsteam(istringstream &istr, string str){
  istr.str(str); // Set the string
  istr.clear(); // Claer states, see https://stackoverflow.com/questions/2767298/c-repeatedly-using-istringstream for the detailed reasons.
}

// Read in the atoms' coordinates
void readposcar_atom_coords(vector<string> &buffer, int &atoms_num, 
                            double* &seq_atom_coords){
  // - Parameters
  istringstream buffer_atom_posi; // Buffer for process Position of the atoms
  int i_line;
  seq_atom_coords = new double[3*atoms_num];
  // - Loop of the atoms coordinates
  for (int i_atom=0; i_atom<atoms_num; ++i_atom){
    i_line = get_xth_atom_line_num_in_poscar(i_atom);
    str_into_istringsteam(buffer_atom_posi, buffer[i_line]);
    // -- Loop for the a1,a2,a3 fractional coordinate
    for (int i_123=0; i_123<3; ++i_123){
      buffer_atom_posi >> seq_atom_coords[3*i_atom + i_123];
    }
  }
}


// Read in the structure data from the POSCAR file
int readin_data_from_poscar(string &file_info, string &scale_factor,
                            double lattice_vec[3][3], int* &type_list, 
                            double* &seq_atom_coords, int &atoms_num,
                            string filename="POSCAR"){
  ifstream poscar_file(filename);
  vector<string> buffer;
  // - Read in all content in the file
  if (poscar_is_not_open(poscar_file, filename)) return 1;
  readlines(poscar_file, buffer); // Read files all content to `buffer`.
  if (poscar_is_not_valid(buffer)) return 1;
  // - Read in the POSCAR info.s
  readposcar_header(buffer, file_info, scale_factor);
  readposcar_lattice(buffer, lattice_vec);
  readposcar_atom_type(buffer, type_list, atoms_num);
  readposcar_atom_coords(buffer, atoms_num, seq_atom_coords);
  // - Program succeed
  return 0;
}


// Use spglib for analysis the atomic strucuture
void spglib_run(double lattice_vec[3][3], double atom_coords[][3],
                int* type_list, int atoms_num, int aperiodic_axis, 
                double symprec, 
                SpglibDataset* &spg_dataset, 
                char ptg_symbol[6], int ptg_trans_mat[3][3]){
  // - Create the spglib database var. and get commom info about the crystal
  spg_dataset = spg_get_layer_dataset(
    lattice_vec, atom_coords, type_list, atoms_num, aperiodic_axis, symprec);
  // - Get point group
  spg_get_pointgroup(
    ptg_symbol, ptg_trans_mat, 
    spg_dataset->rotations, spg_dataset->n_operations);
}

// Decided the space group using spglib
int spglib_2d_tovasp(SpglibDataset* &spg_dataset,
                     char ptg_symbol[6], int ptg_trans_mat[3][3],
                     string filename="POSCAR", double symprec=1E-6){
  // - Parameters
  int _prog_failed;
  string file_info, scale_factor;
  double lattice_vec[3][3];
  int* type_list;
  int atoms_num;
  double* seq_atom_coords; // Atomic coord. 1D array, each 3 number is one coord.
  // - Read in the info.s from POSCAR
  _prog_failed = 
    readin_data_from_poscar(
      file_info, scale_factor, lattice_vec, type_list, 
      seq_atom_coords, atoms_num, filename);
  if (_prog_failed) return 1;
  // - Using spglib
  int aperiodic_axis = 2; // 0,1,2 for x,y,z axis.
  double (*atom_coords)[3] = (double(*)[3])(seq_atom_coords); // Transfer the 1D array to 2D
  spglib_run(
    lattice_vec, atom_coords, type_list, atoms_num, aperiodic_axis, symprec,
    spg_dataset, ptg_symbol, ptg_trans_mat);
  // - Program succeed
  return 0;
}


// Print Help info
void print_help_info(){
  cout << "[help] ./spglib_2d_tovasp.x -f <filename> -p <symprec>" << endl;
}

// Read in the command line parameters
int readin_commandline(int argc, char *argv[], 
                       string &poscar_filename, float &symprec){
  string curr_arg;
  for (int i_arg=1; i_arg<argc; ++i_arg){
    curr_arg = (string)(argv[i_arg]);
    // Determine which option
    if (curr_arg == "-h"){ // help info.
      print_help_info();
      return 1;
    } else if (curr_arg == "-f"){ // filename of the poscar
        if (i_arg < argc-1) poscar_filename = (string)(argv[i_arg+1]);
    } else if (curr_arg == "-p"){ // symmetry precsion
        if (i_arg < argc-1) symprec = stof(argv[i_arg+1]);
    }
  }
  return 0;
}

int main(int argc, char *argv[]){
  int _prog_failed;
  // ***
  // Read in necessary parameters from conmmand line
  // ***
  string poscar_filename = "POSCAR";
  float symprec = 1E-6;
  _prog_failed = readin_commandline(argc, argv, poscar_filename, symprec);
  if (_prog_failed) return 1;
  // ***
  // Analysis the material in POSCAR with spglib
  // ***
  // - Info.s need to be read out
  SpglibDataset* spg_dataset;
  char ptg_symbol[6]; // Point group symbol
  int ptg_trans_mat[3][3]; // Point group transition matrix
  // - Read info.s from POSCAR and analysis it
  _prog_failed = spglib_2d_tovasp(
    spg_dataset, ptg_symbol, ptg_trans_mat,
    poscar_filename, symprec);
  if (_prog_failed) return 1;
  // ***
  // Print result info.
  // ***
  cout << "----------------------------- " << endl;
  cout << "[info] Stru. File  : " << poscar_filename << endl;
  cout << "[info] Symm. Prec. : " << symprec << endl;
  cout << "[info] Space Group : " << spg_dataset->international_symbol << endl;
  cout << "[info] Point Group : " << ptg_symbol << endl;
  cout << "----------------------------- " << endl;
  cout << "[tips] More info. can be print by modify `here` in code..." << endl;
  // ***
  // Free the dataset var.
  // ***
  spg_free_dataset(spg_dataset);
  return 0;
}
