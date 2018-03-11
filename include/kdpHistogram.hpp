#ifndef KDP_HISTOGRAM
#define KDP_HISTOGRAM

#include "kdpTools.hpp"

#include <vector>
#include <set>
#include <string>
#include <stdexcept>

namespace kdp
{
	
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ____                        
    |  _ \ __ _ _ __   __ _  ___ 
    | |_) / _` | '_ \ / _` |/ _ \
    |  _ < (_| | | | | (_| |  __/
    |_| \_\__,_|_| |_|\__, |\___|
                      |___/      
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
		
// A simple class for specifying a range
// A useful feature is the ability to construct from an initalizer_list
// (i.e. Range x = {1, 2})
// When using the ctor, we check that min <= max and
// both min and max are numeric (finite or infinite, but NOT nan)
template<typename edge_t>
struct Range
{
   edge_t min, max;

   Range();
   Range(edge_t const min_in, edge_t const max_in);
   Range(std::initializer_list<edge_t> const& init);

   bool operator == (Range const& that) const;
   bool operator not_eq (Range const& that) { return not(*this == that);}

   private:
      void CheckSortedAndNumeric();
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _____ _       _ _       ____                        
    |  ___(_)_ __ (_) |_ ___|  _ \ __ _ _ __   __ _  ___ 
    | |_  | | '_ \| | __/ _ \ |_) / _` | '_ \ / _` |/ _ \
    |  _| | | | | | | ||  __/  _ < (_| | | | | (_| |  __/
    |_|   |_|_| |_|_|\__\___|_| \_\__,_|_| |_|\__, |\___|
                                              |___/      
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
                                              
// A quick extension to Range, with the ctor requiring min and max be finite
template<typename edge_t>
struct FiniteRange : public Range<edge_t>
{
   FiniteRange():Range<edge_t>() {}   
   FiniteRange(edge_t const min_in, edge_t const max_in);         
   FiniteRange(std::initializer_list<edge_t> const& init);

   private:
      void CheckFinite();
};	

// In this kdpHistogram scheme: 
//  - A "bin" is defined as the semi-inclusive interval [lower, upper).
//  - The correct bin for value x is found by incrementing bins until
//       (x >= lower_i) == false
//    Thus x belongs to the previous bin (i - 1).
// Since only floating point numbers are binned,
// (for which comparison is transitive) this definition is acceptable.

// However, since the values being binned have error
// (either from measurement or floating point rounding error and cancellation),
// there does arise the question of a value which falls
// within Order(error) of an edge -- which bin does it belong to?
// Perhaps its weight should be split between both.
// Since this is generally not a large problem,
// it is currently not addressed in this namespace.

// The floating point type used to define bins (via their edges)
// and therefore also used for the values being binned.
typedef double edge_t;

// The floating point type used to store the weight in each bin
// (since different x will not necessarily have the same weight)
typedef double weight_t;

// Histograms use "bin specs" to define their bins. There are three types:
//    1. UNIFORM ... bins have uniform width:
//       {0, 2, 4, 6, 8, ..., 128}
//    2. LOG_UNIFORM ... bins have uniform width in LOG-space:
//       {0, 1, 2, 4, 8, 16, ..., 128}
//    3. NON_UNIFORM ... bins have non inform width:
//       {0, 32, 48, 56, 60, 62, 64, 66, 68, 72, 80, 96, 128}
enum class BinType {UNIFORM, LOG_UNIFORM, NON_UNIFORM};

// When constructing uniform bins (normal or log),
// the user must decide whether the lowest and highest edge will be
// "pinned" to the supplied range. For example,
// when spanning [0, 1) with bins of width 0.15,
// an integer number of bins will not fit.
// Therefore, if one wants to use a bin width of exactly 0.15,
// either the lowest edge must move below 0 or the highest edge
// must move above 1. The user must decide which edge to "pin".
// Conversely, when spanning [0, 1) with 5 bins, it is obvious
// that BOTH the lowest and highest edges should be "pinned":
//    {0, 0.2, 0.4, 0.6, 0.8, 1.0}
// Returning to spanning [0, 1) with a bin width of exactly 0.15,
// the user could still decide that both edges should be pinned,
// which requires the bin width to be adjusted.
// Currently, bin widths are adjusted by rounding to the nearest
// integer number of bins which span the range.
enum class PinnedEdge {BOTH, LEFT, RIGHT};
// FIX: RIGHT is broken for LOG_UNIFORM

// The BinSpecs class allows multiple histograms to use the
// same bin specs (e.g. the 3 histograms binning the pT of
// the first 3 jets can share the same pT range and bin width).

// UNIFORM and LOG_UNIFORM bins can be defined by the number of bins
// UNIFORM bins can be defined by their width
// LOG_UNIFORM bins can be defined by their "scale"
//    (the ratio of adjacent edges, scale = exp(uniform log Width)).

// BinSpecs creates a list of edges, which ultimately determine
// which bin each value belongs belongs to. To create under/over-flow bins,
// this list is automatically padded with -inf and inf.

// BinSpecs must have a name for a simple reason;
// when the ctor sanity checks fails, and the exception is thrown,
// you'll want to find the offending line in the code (via the name string).

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ____  _       ____                      
    | __ )(_)_ __ / ___| _ __   ___  ___ ___ 
    |  _ \| | '_ \\___ \| '_ \ / _ \/ __/ __|
    | |_) | | | | |___) | |_) |  __/ (__\__ \
    |____/|_|_| |_|____/| .__/ \___|\___|___/
                        |_|                  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
// A set of specs defining a histograms bins
class BinSpecs
{
   private:
      static const edge_t INF; // infinity, for padding with under/over-flow

   public:
      std::string name;
      size_t numBins;  // Does not count the under/over-flow bins
      kdp::FiniteRange<edge_t> range;     // lowest and highest finite edge
      edge_t binWidth; // in linear (UNIFORM) or log space (LOG_UNIFORM)
      BinType binType;     
      PinnedEdge pinned;
      std::vector<edge_t> edges; // List of bin edges

      static const size_t MAX_BINS = 1<<16; // 2^16 = 65,536

      // The ctor for UNIFORM and LOG_UNIFORM bins, given numBins
      BinSpecs(std::string const& name_in,
         size_t const numBins_in,
         std::initializer_list<edge_t> const& range_init,
         BinType const binType_in = BinType::UNIFORM);

      // The ctor for UNIFORM and LOG_UNIFORM bins, given width/scale (respectively).
      // Note that the width/scale occurs AFTER the range,
      // to unambiguously differentiate from the numBins ctor.
      BinSpecs(std::string const& name_in,
         std::initializer_list<edge_t> const& range_init,
         edge_t widthOrScale,
         BinType const binType_in = BinType::UNIFORM,
         PinnedEdge const pinned_in = PinnedEdge::LEFT);

      // The (only) ctor for NON_UNIFORM bins
      BinSpecs(std::string const& name_in, std::vector<edge_t> const& edges_in);
      
      ~BinSpecs();

      bool operator == (BinSpecs const& that) const;
      bool operator != (BinSpecs const& that) const 
         {return not (this->operator==(that));}

   private:
      void SharedSanityChecks(); // Used by every ctor
      void UniformSanityChecks(); // Used by UNIFORM and LOG_UNIFORM
      void CalculatedNumBinsSanityChecks(); // Used after Fill_Uniform_Edges
      void FillEdges(bool const fromNumBins);
      void Fill_Uniform_Edges(bool const fromNumBins);
      void ExceptionHandler(std::exception& except);
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    __        __   _       _     _   _____                     
    \ \      / /__(_) __ _| |__ | |_| ____|_ __ _ __ ___  _ __ 
     \ \ /\ / / _ \ |/ _` | '_ \| __|  _| | '__| '__/ _ \| '__|
      \ V  V /  __/ | (_| | | | | |_| |___| |  | | | (_) | |   
       \_/\_/ \___|_|\__, |_| |_|\__|_____|_|  |_|  \___/|_|   
                     |___/                                     
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

// Placing weight and error in the same structure speeds up binning,
// because both values will reside on the same cache line
// (as opposed to separate vectors storing weight and error2)
struct WeightError
{
   weight_t weight;
   weight_t error2;

   WeightError():weight(0.), error2(0.) {}
   WeightError(weight_t const weight_in, weight_t const error2_in):
      weight(weight_in), error2(error2_in) {}

   // We can divide WeightError (mostly useful for efficiency plots)
   WeightError& operator/=(WeightError const& that);
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _     __     __        
    | |/ /  _ \|  _ \| | | (_)___| |_ __\ \   / /__  ___ 
    | ' /| | | | |_) | |_| | / __| __/ _ \ \ / / _ \/ __|
    | . \| |_| |  __/|  _  | \__ \ || (_) \ V /  __/ (__ 
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \_/ \___|\___|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

// A generic vector to store binned data.
// This makes it easier to create and manipulate multi-dimensional histograms
template<class content>
class KDPHistoVec : public std::vector<content>
{
   public:
      // Each KDPHistoVec is a std::vector (simplifying manipulation)
      // It must contain a reference to some BinSpecs
      BinSpecs const& binSpecs;
      // It contains an iterator to the edge which stores the min
      // (this elludes a branch every time FindBin is called)
      decltype(binSpecs.edges.begin()) const minEdge;    

   private:
      // Single final ctor, so the minEdge code is only in one place
      // Steals vector, links to binSpecs, redefines minEdge (just in case)
      KDPHistoVec(std::vector<content>&& disowned,
         BinSpecs const& binSpecs_in):
      std::vector<content>(std::move(disowned)),
         binSpecs(binSpecs_in),
         minEdge(binSpecs.edges.begin() +
            ((binSpecs.binType == BinType::LOG_UNIFORM) ? 2 : 1)) {}
      
   public:
      // Normal ctor; empty KDPHistoVec is created recursively,
      // given a list of BinSpecs (one for each dimension)
      template<typename ...Args>
      KDPHistoVec(BinSpecs const& binSpecs_in, Args && ...otherBinSpecs):
      // Always add 2 bins for under and overflow
      KDPHistoVec(std::move(std::vector<content>(binSpecs_in.numBins + 2,
         content(std::forward<Args>(otherBinSpecs)...))), binSpecs_in) {}

      // Copy ctor copies incoming vector
      KDPHistoVec(KDPHistoVec<content> const& other):
         KDPHistoVec(std::move(std::vector<content>(other)), other.binSpecs) {}

      // Move constructor steals incoming vector
      KDPHistoVec(KDPHistoVec<content> && disowned):
         KDPHistoVec(std::move(disowned), disowned.binSpecs) {}

      // Given a value, find the "bin" it belongs to
      // (in a 1D array, content = WeightError)
      // (in a 2D array, content = KDPHistoVec<WeightError>)
      content& FindBin(edge_t const val);

      // We can divide two histovecs by each other (for efficiency histograms)
      KDPHistoVec& operator/=(KDPHistoVec<content> const& that);
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _   _                            _ _         
    | \ | | ___  _ __ _ __ ___   __ _| (_)_______ 
    |  \| |/ _ \| '__| '_ ` _ \ / _` | | |_  / _ \
    | |\  | (_) | |  | | | | | | (_| | | |/ /  __/
    |_| \_|\___/|_|  |_| |_| |_|\__,_|_|_/___\___|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/      

// Normalize is an enum to handle different histogram normalization schemes
//    NO: don't perform any normalization, output the binned weight verbatim
//    INTEGRATABLE: divide bin weight by bin width
//    UNITARY: divide bin weight by total binned weight (including under/over-flow)
//    PDF: both INTEGRATABLE and UNITARY (i.e. like a probability distribution)
//       NOTE: PDF = INTEGRATABLE & UNITARY
//       (i.e. it's not 3 because 3 comes after 2)

typedef uint_fast8_t norm_t;

enum class Normalize : norm_t {NO = 0, INTEGRATABLE = 1, UNITARY = 2, PDF = 3};

// In order to use the type-safe enum class Normalize like a bit flag,
// we have to  manually define binary operations.
// Must make sure their declaration as kdpHistogram::operator,
// or you will get undefined reference errors during linking
Normalize operator & (Normalize const lhs, Normalize const rhs);
Normalize operator | (Normalize const lhs, Normalize const rhs);
   
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _                                      _                    
    | |/ /  _ \|  _ \| | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___     | |__   __ _ ___  ___ 
    | ' /| | | | |_) | |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \    | '_ \ / _` / __|/ _ \
    | . \| |_| |  __/|  _  | \__ \ || (_) | (_| | | | (_| | | | | | |   | |_) | (_| \__ \  __/
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_|___|_.__/ \__,_|___/\___|
                                           |___/                   |_____|                    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

// The base class for multi-dimensional histograms, with basic housekeeping,
// but without the actual histogram (avoid too much templating,
// because behavior changes too much between 1D and 2D implementations)
// It ensures that the fileName is unique, so only one active histogram
// (per thread, not thread safe) is writing to a given file
// Inheritied classes only have to define Fill functions and Write
class KDPHistogram_base
{
   private:
       // A static list of active names, to ensure files are unique
      static std::set<std::string> fileNameMutex;
      
   protected:
      std::string fileName;
      
      weight_t totalWeight; // Total weight binned so far
      weight_t writeScale; // Multiply by all weights when writing
      Normalize norm; 
      bool written;

      // protected ctor; called by inheriting classes
      KDPHistogram_base(std::string const& fileName_in,
         weight_t const totalWeight_in,
         weight_t const writeScale_in,
         Normalize const norm_in);

      // Machinery to throw exceptions identifying the histo by fileName
      std::string ExceptionPrefix() const;
      void RTExcept(std::string const& info) const;

      // Make sure weight and error aren't nan or inf
      inline void CheckWeightError(weight_t const weight, weight_t const error) const;

      // Called by the destructors of derived classes
      void AutomaticWrite();

      // Write the common header to the histo file
      void WriteHeader(std::ofstream& outFile, std::string const& className);
   public:
      void SetWriteScale(weight_t const newScale) {writeScale = newScale;}
      void SetNormalize(Normalize const newNorm) {norm = newNorm;}
      
      // Write the histogram (depends on dimensionality)
      virtual void Write() = 0;
      void Write(weight_t const scale) {SetWriteScale(scale); this->Write();}         
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _                                  _ 
    | |/ /  _ \|  _ \| | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___ / |
    | ' /| | | | |_) | |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \| |
    | . \| |_| |  __/|  _  | \__ \ || (_) | (_| | | | (_| | | | | | | |
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_|_|
                                           |___/                       
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

// 1D histogram
//    gnuplot> plot "my.histo1D" using 1:2:3 with errorbars
class KDPHistogram1 : public KDPHistogram_base
{
   private:
      KDPHistoVec<WeightError> hist;
      static const std::string FILE_EXTENSION;

      // One ctor, called by all others
      KDPHistogram1(std::string const& name,
         weight_t const totalWeight_in,
         weight_t const writeScale_in,
         Normalize const norm_in,
         KDPHistoVec<WeightError>&& disowned);

   public:
      // ctor given writeScale and norm
      KDPHistogram1(std::string const& name,
         BinSpecs const& xBinSpecs,
         weight_t const writeScale_in,
         Normalize const norm_in = Normalize::NO);  
   
      // Minimal ctor (writeScale = 1., norm = PDF)
      KDPHistogram1(std::string const& name,
         BinSpecs const& xBinSpecs,
         Normalize const norm_in = Normalize::PDF);            

      ~KDPHistogram1() {AutomaticWrite();}

      void Fill(edge_t const x, weight_t const weight = 1.) {Fill(x, weight, weight);}
      void Fill(edge_t const x, weight_t const weight, weight_t const error);

      // Divide two histograms to create (and immedietely write)
      // a histogram storing their ratio
      void WriteRatio(KDPHistogram1 const& that, std::string const& newName) const;
      
      void Write();         
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _                                  ____  
    | |/ /  _ \|  _ \| | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___ |___ \ 
    | ' /| | | | |_) | |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \  __) |
    | . \| |_| |  __/|  _  | \__ \ || (_) | (_| | | | (_| | | | | | |/ __/ 
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_|_____|
                                           |___/                           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

// A 2D histogram
//    gnuplot> set pm3d corners2color c1
//    gnuplot> splot "my.histo2D" with pm3d
class KDPHistogram2 : public KDPHistogram_base
{
   private:
      KDPHistoVec<KDPHistoVec<WeightError> > hist;
      static const std::string FILE_EXTENSION;

      KDPHistogram2(std::string const& name,
         weight_t const totalWeight_in,
         weight_t const writeScale_in,
         Normalize const norm_in,
         KDPHistoVec<KDPHistoVec<WeightError> >&& disowned);

         void WriteUnderOverX(std::ofstream& outFile, weight_t const normWeight);
         void WriteUnderOverY(std::ofstream& outFile, weight_t const normWeight);

   public:
      // ctor given writeScale and norm
      KDPHistogram2(std::string const& name,
         BinSpecs const& xBinSpecs, BinSpecs const& yBinSpecs,
         weight_t const writeScale_in,
         Normalize const norm_in = Normalize::NO);  
   
      // Minimal ctor (writeScale = 1., norm = PDF)
      KDPHistogram2(std::string const& name,
         BinSpecs const& xBinSpecs, BinSpecs const& yBinSpecs,
         Normalize const norm_in = Normalize::PDF);            
      
      ~KDPHistogram2() {AutomaticWrite();}

      // Fill with weight and weight error
      void Fill(edge_t const x, edge_t const y, weight_t const weight = 1.)
         {Fill(x, y, weight, weight);}
         
      void Fill(edge_t const x, edge_t const y,
         weight_t const weight, weight_t const error);

      // Divide two histograms to create (and immedietely write)
      // a histogram storing their ratio
      void WriteRatio(KDPHistogram2 const& that, std::string const& newName) const;

      void Write();
};

}

/*
// A simple file mutex, to prevent multiple objects/threads
// from concurrent access to a file.
// Each method returns 0 if the operation was successful,
// and a non-zero number otherwise.
class FileMutex
{
   private:
      static std::mutex fileAccess;

   public:
      static const std::string EXTENSION;

      static int Lock(std::string const& fileName);
      static int Unlock(std::string const& fileName);         
};

std::string const FileMutex::EXTENSION = ".mutex";
std::mutex FileMutex::fileAccess;

int FileMutex::Lock(std::string const& fileName)
{
   std::string const mutexFileName = fileName + EXTENSION;
   std::fstream fileMutex;

   // Make sure this thread is the only with file access by locking fileAccess,
   // then automatically unlock fileAccess upon function return.
   std::lock_guard<std::mutex> safeLock(fileAccess);

   // Try to open file for input; a non-existent file will not be created
   fileMutex.open(mutexFileName, std::ios::in);

   // If the file is open, then it already exists. Return error
   if(fileMutex.is_open()) return 1;

   // Otherwise the file doesn't exist, create it to block other KDPHistograms
   fileMutex.open(mutexFileName, std::ios::out);
   return 0;
}

int FileMutex::Unlock(std::string const& fileName)
{
   std::lock_guard<std::mutex> safeLock(fileAccess);
   return std::remove((fileName + EXTENSION).c_str());
}

*/


#endif
