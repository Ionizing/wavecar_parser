/*
 * Plot KS pseudo wavefunction at specific kpoint, specific band.
 * An example that uses wavehigh class
 */

#include <wavehigh.hpp>
#include <cxxopts.hpp>

using ionizing::WAVEHIGH;

int main(int argc, char* argv[]) {
  // int argc_bak = argc;
  cxxopts::Options options{argv[0],"A wavefunction plotter using WAVECAR."};

  bool 
    is_verbose       = false,
    is_real          = false,
    is_print_help    = false,
    is_list_info     = false,
    is_print_example = false;

  int  
    ikpoint    = 0,
    iband      = 0,
    ispin      = 0;

  std::string prefix{"wfc"};
  std::string wavecar_fname{"WAVECAR"};
  std::string poscar_fname{"POSCAR"};

  options
    .allow_unrecognised_options()
    .add_options()
    // ("v, verbose", "Enable verbose output for debugging .",
        // cxxopts::value<bool>(is_verbose))
    ("b, band"   , "Specify which band to plot,            default: 0",
        cxxopts::value<int>(iband))
    ("e, example", "Give example usage of plotWave",
        cxxopts::value<bool>(is_print_example))
    ("h, help", "Print help",
        cxxopts::value<bool>(is_print_help))
    ("k, kpoint" , "Specify which kpoint to plot,          default: 0",
        cxxopts::value<int>(ikpoint))
    ("l, list"   , "List brief info of current wavecar",
        cxxopts::value<bool>(is_list_info))
    ("prefix"    , "File name prefix of output .vasp file, default: \"wfc\"",
        cxxopts::value<std::string>(prefix))
    ("poscar"    , "Specify which POSCAR to read,      default: \"POSCAR\"",
        cxxopts::value<std::string>(poscar_fname))
    ("r, real"   , "Only output real part of wavefunction",
        cxxopts::value<bool>(is_real))
    ("s, spin"   , "Specify which spin to plot,            default: 0",
        cxxopts::value<int>(ispin))
    ("w, wavecar", "Specify which WAVECAR to parse,    default: \"WAVECAR\"", 
        cxxopts::value<std::string>(wavecar_fname))
    ;

  auto result = options.parse(argc, argv);
  is_verbose = false; // unused, will impl 100 yrs later

  if (is_print_help /* or 1 == argc_bak */) {
    std::cout << options.help();
    return 0;
  }

#if defined DEBUG_ARGPARSE
  { // for debug
    std::cout << std::boolalpha 
      << "is_real = " << is_real << std::endl
      << "ispin = " << ispin << std::endl
      << "ikpoint = " << ikpoint << std::endl
      << "iband = " << iband << std::endl
      << "prefix = " << prefix << std::endl
      << "wavecar = " << wavecar_fname << std::endl
      << "is_list_info = " << is_list_info << std::endl
      << "is_print_help = " << is_print_help << std::endl
      << "is_print_example = " << is_print_example << std::endl
      ;
  }
#endif

  if (is_print_example) {
    std::cout << "Example:" << std::endl
     << "    For wave plot: " << argv[0] 
        << " --spin=0 --kpoint=0 --band=0 --poscar=POSCAR --prefix=wfc --wavecar=WAVECAR --real" << std::endl
     << "       Or shortly: " << argv[0] 
        << " -s 0 -k 0 -b 0 --poscar POSCAR --prefix wfc -w WAVECAR -r" << std::endl
     << "    For wave info: " << argv[0] << " --list or " << argv[0] << " -l" << std::endl
     ;

    return 0;
  }

  WAVEHIGH wave{wavecar_fname.c_str()}; 

  if (is_list_info) {
    wave.printInfo(std::cout);
  } else {
    wave.plotWave(ispin, ikpoint, iband, poscar_fname.c_str(), prefix.c_str(), is_real);
  }
  return 0;
}
