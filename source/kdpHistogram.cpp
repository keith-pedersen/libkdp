// Simple histogrammer

#include "kdpHistogram.hpp"

#include <fstream>
#include <limits>
#include <iostream>
#include <assert.h>

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ____                        
    |  _ \ __ _ _ __   __ _  ___ 
    | |_) / _` | '_ \ / _` |/ _ \
    |  _ < (_| | | | | (_| |  __/
    |_| \_\__,_|_| |_|\__, |\___|
                      |___/      
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

template<typename edge_t>
kdp::Range<edge_t>::Range(): min(0.), max(1.) {}

template<typename edge_t>
kdp::Range<edge_t>::Range(edge_t const min_in, edge_t const max_in):
   min(min_in), max(max_in)
{CheckSortedAndNumeric();}

template<typename edge_t>
kdp::Range<edge_t>::Range(std::initializer_list<edge_t> const& init)
{
   if(init.size() not_eq 2)
      throw std::invalid_argument("|Range| initializer_list must have 2 values");

   min = *init.begin();
   max = *(init.begin() + 1);

   CheckSortedAndNumeric();
}

template<typename edge_t>
bool kdp::Range<edge_t>::operator == (Range const& that) const
{
   return (this->min == that.min) and (this->max == that.max);
}

template<typename edge_t>
void kdp::Range<edge_t>::CheckSortedAndNumeric()
{
   if(min > max) throw std::invalid_argument("|Range| min > max");
   if(std::isnan(min)) throw std::invalid_argument("|Range| min = nan");
   if(std::isnan(max)) throw std::invalid_argument("|Range| max = nan");
}

template class kdp::Range<float>;
template class kdp::Range<double>;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _____ _       _ _       ____                        
    |  ___(_)_ __ (_) |_ ___|  _ \ __ _ _ __   __ _  ___ 
    | |_  | | '_ \| | __/ _ \ |_) / _` | '_ \ / _` |/ _ \
    |  _| | | | | | | ||  __/  _ < (_| | | | | (_| |  __/
    |_|   |_|_| |_|_|\__\___|_| \_\__,_|_| |_|\__, |\___|
                                              |___/      
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
template<typename edge_t>
kdp::FiniteRange<edge_t>::FiniteRange(edge_t const min_in, edge_t const max_in):
   Range<edge_t>(min_in, max_in)
{CheckFinite();}

template<typename edge_t>
kdp::FiniteRange<edge_t>::FiniteRange(std::initializer_list<edge_t> const& init):
   Range<edge_t>(init)
{CheckFinite();}

template<typename edge_t>
void kdp::FiniteRange<edge_t>::CheckFinite()
{
   if(std::isinf(this->min)) throw std::invalid_argument("|Range| min = -inf, must be finite");
   if(std::isinf(this->max)) throw std::invalid_argument("|Range| max = inf, must be finite");
}

template class kdp::FiniteRange<float>;
template class kdp::FiniteRange<double>;



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ____  _       ____                      
    | __ )(_)_ __ / ___| _ __   ___  ___ ___ 
    |  _ \| | '_ \\___ \| '_ \ / _ \/ __/ __|
    | |_) | | | | |___) | |_) |  __/ (__\__ \
    |____/|_|_| |_|____/| .__/ \___|\___|___/
                        |_|                  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

const kdp::edge_t kdp::BinSpecs::INF =
   std::numeric_limits<edge_t>::infinity();

// numBins ctor
kdp::BinSpecs::BinSpecs(std::string const& name_in,
   size_t const numBins_in,
   std::initializer_list<edge_t> const& range_init,
   BinType const binType_in):
name(name_in), numBins(numBins_in),
range(), // Set below, to catch and rethrow exceptions, appending name
binWidth(-1.), // Set below, to keep conditional out of init list
binType(binType_in),
pinned(PinnedEdge::BOTH) // Only option, given numBins
{
   try
   {
      if(numBins == 0) throw std::invalid_argument("not enough bins supplied");
      range = range_init;

      UniformSanityChecks();

      binWidth = ((binType == BinType::UNIFORM) ?
      (range.max - range.min) : std::log(range.max/range.min))/edge_t(numBins);

      FillEdges(false);
   }
   catch(std::exception& e)
   {
      ExceptionHandler(e);
   }
}

// width/scale ctor
kdp::BinSpecs::BinSpecs(std::string const& name_in,
   std::initializer_list<edge_t> const& range_init,
   edge_t widthOrScale,
   BinType const binType_in, PinnedEdge const pinned_in):
name(name_in),
numBins(1), // Default value, will determine in FillEdges
range(), // Set below, to catch and rethrow exceptions, appending name
binWidth(-1.), // Set below, to keep conditional out of init list
binType(binType_in), pinned(pinned_in)
{
   try
   {
      range = range_init;
      UniformSanityChecks();

      binWidth = (binType_in == BinType::UNIFORM) ?
         widthOrScale : std::log(widthOrScale);
      
      FillEdges(true);
   }
   catch(std::exception& e)
   {
      ExceptionHandler(e);
   }
}   

// NON_UNIFORM ctor
kdp::BinSpecs::BinSpecs(std::string const& name_in,
   std::vector<edge_t> const& edges_in):
   name(name_in),
   numBins(edges_in.size() - 1), // this could be very large if edges is empty
   range(),// If edges_in is empty, trying to construct range would cause SEG_FAULT
   binWidth(-1.), // nonsensical value
   binType(BinType::NON_UNIFORM),
   pinned(PinnedEdge::BOTH) // Only option, given edges
{
   try
   {
      if(edges_in.size() < 2) throw std::invalid_argument("< 2 edges supplied, cannot form a bin");
      else range = {edges_in.front(), edges_in.back()};

      SharedSanityChecks();

      // Verify the supplied edges are all finite
      // DO NOT attempt to fix anything, let the user decide how to fix it
      // (that way they are not surprised about this behavior under the hood)
      {
         auto prevEdge = edges_in.begin(); 
         for(auto thisEdge = prevEdge + 1; thisEdge not_eq edges_in.end(); ++prevEdge, ++thisEdge)
         {
            // The duplicate check will only find all duplicates in a sorted list
            // But an unsorted list will fail the sorting test
            if(*thisEdge == *prevEdge) throw std::invalid_argument("edge list has duplicates");
            if(*thisEdge < *prevEdge) throw std::invalid_argument("edge list not sorted");         
            if(std::isinf(*thisEdge)) throw std::invalid_argument("supplied edge is +/-inf");
            if(std::isnan(*thisEdge)) throw std::invalid_argument("supplied edge is nan");         
         }
      }

      // Copy in input edges, padding infinity on each side
      edges.reserve(numBins + 2);
      edges.push_back(-INF);
      edges.insert(edges.end(), edges_in.begin(), edges_in.end());
      edges.push_back(INF);
   }
   catch(std::exception& e)
   {
      ExceptionHandler(e);
   }
}

// Define this, otherwise -Winline complains about not inlining it
kdp::BinSpecs::~BinSpecs(){}

// Sanity checks used by all three ctors
void kdp::BinSpecs::SharedSanityChecks()
{
   if(numBins > BinSpecs::MAX_BINS)
      throw std::invalid_argument(std::to_string(numBins) + " = too many bins supplied");
}

// Sanity checks only used by UNIFORM/NON_UNIFORM
void kdp::BinSpecs::UniformSanityChecks()
{
   if(binType == BinType::NON_UNIFORM)
   {
      // wrong ctor used
      throw std::invalid_argument("for NON_UNIFORM bins, edges must be supplied in a vector");
   }
   else if(range.min == range.max) // FiniteRange will accept this, but we can't
      throw std::invalid_argument("min == max, cannot create bins");
      
   SharedSanityChecks();
}

// When width/scale is supplied, make sure the resulting numBins is sensible
void kdp::BinSpecs::CalculatedNumBinsSanityChecks()
{
   if(numBins == 0)
      throw std::invalid_argument("Zero bins from width/scale ... this shouldn't happen (code logic error)");
   else if(numBins > BinSpecs::MAX_BINS)
      throw std::invalid_argument(std::to_string(numBins) + " = too many bins from width/scale");
}

bool kdp::BinSpecs::operator==(BinSpecs const& that) const
{
   if(this == &that) return true;
   else
   {
      // There must be the same number of bins spanning the same range
      if((this->numBins == that.numBins) and
         (this->range == that.range))
      {
         // If either is NON_UNIFORM, the only way to be absolutely sure
         // is to check edge by edge
         if((this->binType == BinType::NON_UNIFORM) or (that.binType == BinType::NON_UNIFORM))
         {
            auto itThis = this->edges.begin();
            auto itThat = that.edges.begin();

            // Keep iterating while both edges are the same and we
            // haven't hit the end of either list
            while(
               (itThis not_eq this->edges.end()) and
               (itThat not_eq that.edges.end()) and
               (*itThis == *itThat))
            {++itThis; ++itThat;}

            // If all edges are the same, then both iterators should be
            // at the end of their lists
            return (itThis == this->edges.end()) and (itThat == that.edges.end());
         }
         else
         {
            // Otherwise, if both are uniform, then as long as they have
            // identical binWidths and binTypes, they're the same
            // (even if pinned is different). This is due to the algorithm
            // which fills the bins, since it only depends on
            // min, max, numBins, binWidth and binType
            return ((this->binWidth == that.binWidth) and 
               (this->binType == that.binType));
         }
      }
      else return false;
   }
}

void kdp::BinSpecs::FillEdges(bool const fromWidth)
{
   switch(binType)
   {
      case BinType::UNIFORM:
         if(binWidth <= 0.) throw std::invalid_argument("bin width must be > 0");

         Fill_Uniform_Edges(fromWidth);
      break;
      
      case BinType::LOG_UNIFORM:
         if(range.min <= 0.)
            throw std::invalid_argument("for LOG_UNIFORM, lower edge must be >= 0");
         if(binWidth <= 0.)
            throw std::invalid_argument("bin scale must be > 1.");

         Fill_Uniform_Edges(fromWidth);
      break;

      default:
         throw std::runtime_error("BinSpecs: no code to handle a new BinType (my mistake)");
      break;
   }
}

// ALWAYS FILL UNIFORM EDGES USING (min*i + max*(N-i)) scheme
// Introduces rounding error in every edge, but it's always Order(epsilon)
// (instead of accumulating as edges are defined incrementally)
void kdp::BinSpecs::Fill_Uniform_Edges(bool const fromWidth)
{
   // Critical ASSUMPTION throughout:
   //    if binType != UNIFORM, then binType == LOG_UNIFORM
   if(fromWidth)
   {
      edge_t& min = range.min;
      edge_t& max = range.max;

      const edge_t rangeWidth = (binType == BinType::UNIFORM) ?
         (max - min) : std::log(max/min);

      if(pinned == PinnedEdge::BOTH)
      {
         // If pinning of both edges requested,
         // change the width by first rounding to an integer number of bins
         numBins = std::max(size_t(std::round(rangeWidth/binWidth)), size_t(1));
         binWidth = rangeWidth/edge_t(numBins);
      }
      else
      {
         // Otherwise, use just enough bins to accommodate rangeWidth
         numBins = std::ceil(rangeWidth/binWidth);

         // Now move min or max 
         if(pinned == PinnedEdge::LEFT)
         {
            max = (binType == BinType::UNIFORM) ?
               (min + numBins*binWidth) : (min * std::exp(numBins*binWidth));
         }
         else if(pinned == PinnedEdge::RIGHT)
         {
            min = (binType == BinType::UNIFORM) ?
               (max - numBins*binWidth) : (max * std::exp(-numBins*binWidth));
         }
      }

      CalculatedNumBinsSanityChecks();
   }

   const edge_t minLoop = (binType == BinType::UNIFORM) ?
      range.min : std::log(range.min);
   const edge_t maxLoop = (binType == BinType::UNIFORM) ?
      range.max : std::log(range.max);

   // Reserve 2 extra slots for -inf/inf
   edges.reserve(numBins + 2);
   // Pad with infinity before/after
   edges.push_back(-INF);

   // For LOG_UNIFORM, add the edge for the zero bin (pseudo-underflow) 
   if(binType == BinType::LOG_UNIFORM) edges.push_back(0.);

   // Manually push_back min/max, to avoid rounding error
   edges.push_back(range.min);

   // Even if using bin width, define edges in terms of num bins
   // This prevents rounding error from accumulating, so even though
   // the bins may have slightly uneven width, they are more likely
   // to exist in the right place.
   {
      size_t minCounter = numBins - 1;
      for(size_t maxCounter = 1; maxCounter < numBins; ++maxCounter, --minCounter)
      {
         edge_t thisEdge = (minLoop*minCounter + maxLoop*maxCounter)/edge_t(numBins);
         if(binType == BinType::LOG_UNIFORM) thisEdge = std::exp(thisEdge);

         edges.push_back(thisEdge);
      }
   }

   edges.push_back(range.max);
   edges.push_back(INF);

   // Account for automatic addition of the zero bin for LOG_UNIFORM
   // (this must be done here, after the loop,
   // otherwise too many edges will be added in the loop)
   if(binType == BinType::LOG_UNIFORM) ++numBins;
}

void kdp::BinSpecs::ExceptionHandler(std::exception& except)
{
   throw std::invalid_argument("BinSpecs: (" + name + ") ... " +
      except.what());
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    __        __   _       _     _   _____                     
    \ \      / /__(_) __ _| |__ | |_| ____|_ __ _ __ ___  _ __ 
     \ \ /\ / / _ \ |/ _` | '_ \| __|  _| | '__| '__/ _ \| '__|
      \ V  V /  __/ | (_| | | | | |_| |___| |  | | | (_) | |   
       \_/\_/ \___|_|\__, |_| |_|\__|_____|_|  |_|  \___/|_|   
                     |___/                                     
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

kdp::WeightError& kdp::WeightError::operator /= (WeightError const& that)
{
   // We need to avoid 0 / 0 ... this way also catches 0/something
   if(this->weight > 0.) 
   {
      // For * and /, *relative* errors add in quadrature
      // distribute for dz^2:   dZ/z = sqrt((dx/x)^2 + (dy/y)^2)
      // Note: a^2/b^2 is more accurate than (a/b)^2
      this->error2 = 
         std::fma(that.error2, kdp::Squared(this->weight)/kdp::Squared(that.weight), this->error2)/
            kdp::Squared(that.weight);
      this->weight = this->weight / that.weight;
   }
   return *this;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _     __     __        
    | |/ /  _ \|  _ \| | | (_)___| |_ __\ \   / /__  ___ 
    | ' /| | | | |_) | |_| | / __| __/ _ \ \ / / _ \/ __|
    | . \| |_| |  __/|  _  | \__ \ || (_) \ V /  __/ (__ 
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \_/ \___|\___|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

template<class content>
content& kdp::KDPHistoVec<content>::FindBin(edge_t const val)
{
   if(val >= binSpecs.range.max) return this->back();
   else
   {
      auto bin = this->begin();

      if(val >= binSpecs.range.min)
      {
         auto leftEdge = minEdge;

         // For UNIFORM/LOG_UNIFORM bins, no need to search every edge
         // However, there might be rounding error while calculating the bin. 
         // Start with the calculated/expected bin and make sure we
         // search past it. This will catch bin shifts to the left or right.

         if(binSpecs.binType == BinType::UNIFORM)
         {
            leftEdge += size_t((val - binSpecs.range.min)/binSpecs.binWidth);
         }
         else if(binSpecs.binType == BinType::LOG_UNIFORM)
         {
            ++bin; // minEdge skipped the zero bin

            // This is actually slower, if there are less than about 100 bins
            // But it has constant execution time, for larger number of bins
            leftEdge += size_t(std::log(val/binSpecs.range.min)/binSpecs.binWidth);            
         }

         // Keep moving the left edge until it is to our right
         // Because +inf is the final value in the edge list,
         // and we already checked for overflow,
         // we do not need an extra check to prevent running off the list

         while(val >= *leftEdge) ++leftEdge;
         bin += size_t(leftEdge - minEdge);
      }
      else if((binSpecs.binType == BinType::LOG_UNIFORM) and (val >= 0.)) ++bin;
      // For LOG_UNIFORM, use the zero bin [0, min)

   return *bin;
   }
}

template<class content>
kdp::KDPHistoVec<content>&
kdp::KDPHistoVec<content>::operator/=(KDPHistoVec<content> const& that)
{
   if(this->binSpecs not_eq that.binSpecs)
      throw std::runtime_error("KDPHistoVec: Cannot divide, incompatible BinSpecs");

   auto itDenom = that.begin();

   for(auto itNumer = this->begin(); itNumer not_eq this->end(); ++itNumer, ++itDenom)
   {
      // This will call this function recursively until content = WeightError
      // This makes it extremely simple to divide multi-D histograms
      (*itNumer)/=(*itDenom);
   }

   return *this;
}

template class kdp::KDPHistoVec<kdp::WeightError>;
template class kdp::KDPHistoVec<kdp::KDPHistoVec<kdp::WeightError> >;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _   _                            _ _         
    | \ | | ___  _ __ _ __ ___   __ _| (_)_______ 
    |  \| |/ _ \| '__| '_ ` _ \ / _` | | |_  / _ \
    | |\  | (_) | |  | | | | | | (_| | | |/ /  __/
    |_| \_|\___/|_|  |_| |_| |_|\__,_|_|_/___\___|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

// Have to declare as kdp:operator, otherwise defined in global scope
kdp::Normalize kdp::operator&(Normalize const lhs,
   Normalize const rhs)
{
   return static_cast<Normalize>
      (static_cast<norm_t>(lhs) & static_cast<norm_t>(rhs));
}

kdp::Normalize kdp::operator|(Normalize const lhs, Normalize const rhs)
{
   return static_cast<Normalize>
      (static_cast<norm_t>(lhs) | static_cast<norm_t>(rhs));
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _                                      _                    
    | |/ /  _ \|  _ \| | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___     | |__   __ _ ___  ___ 
    | ' /| | | | |_) | |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \    | '_ \ / _` / __|/ _ \
    | . \| |_| |  __/|  _  | \__ \ || (_) | (_| | | | (_| | | | | | |   | |_) | (_| \__ \  __/
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_|___|_.__/ \__,_|___/\___|
                                           |___/                   |_____|                    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

std::set<std::string> kdp::KDPHistogram_base::fileNameMutex;

kdp::KDPHistogram_base::KDPHistogram_base(std::string const& fileName_in,
   weight_t const totalWeight_in,
   weight_t const writeScale_in,
   Normalize const norm_in):
fileName(fileName_in),
totalWeight(totalWeight_in), writeScale(writeScale_in), norm(norm_in),
written(false)
{
   // nameMutex resturns false in the 2nd position if the name has already been used
   if(not (fileNameMutex.insert(fileName).second))
      RTExcept("another KDPHistogram is already using this fileName");

   // Clear the file, if it already exists, and ensure we have write privileges
   if(not(std::ofstream(fileName.c_str(), std::ios::out | std::ios::trunc)))
      RTExcept("cannot open file (" + fileName  +") for writing");
   // Now remove the empty file, in case the histogram isn't used
   // (because then it won't be written, which will be easier to detect
   // if there is no file there
   std::remove(fileName.c_str());
}

void kdp::KDPHistogram_base::RTExcept(std::string const& info) const
{
   throw std::runtime_error(ExceptionPrefix() + info);
}

std::string kdp::KDPHistogram_base::ExceptionPrefix() const
{
   return "KDPHistogram: (" + fileName + ") ... ";
}

void kdp::KDPHistogram_base::CheckWeightError(weight_t const weight,
   weight_t const error) const
{
   if(std::isinf(weight)) RTExcept("binning weight is inf");
   if(std::isnan(weight)) RTExcept("binning weight is nan");
   if(weight <= 0.) RTExcept("binning weight is <= 0");

   if(std::isinf(error)) RTExcept("binning error is inf");
   if(std::isnan(error)) RTExcept("binning error is nan");
   if(error < 0.) RTExcept("binning error is < 0");
}

// This must be called from the derived classes, because Write is virtual
void kdp::KDPHistogram_base::AutomaticWrite()
{
   try
   {
      // Don't erase fileName, so we don't screw ourselves within a single run
         // fileNameMutex.erase(fileName);

      // Write if no write was manual requested, and if weight has been binned
      if((not written) and (totalWeight > 0.))
         this->Write(); // This should not throw an exception, but it might
   }
   catch(std::exception const& e)
   {
      std::cerr << "Error during KDPHistogram1 deconstruction (probably Write())\n\n";
      std::cerr << e.what();
   }
}

void kdp::KDPHistogram_base::WriteHeader(std::ofstream& outFile,
   std::string const& className)
{
   const std::string blankLine = "#\n";
   char outBuff[1024]; // should be long enough for one line (ASSUMPTION)
   
   outFile << "# Histogram generated by (" + className + "), ";
   outFile << "in a format easily recognized by gnuplot.\n";
   outFile << blankLine;

   sprintf(outBuff, "# Total weight binned:      %.16e \n", totalWeight);
   outFile << outBuff;
   outFile << blankLine;

   if(bool(norm & Normalize::UNITARY))
   {
      outFile << "# histogram is UNITARY       (bin weights, including under/over, ";
      outFile << "are divided by total weight)\n";
   }

   if(bool(norm & Normalize::INTEGRATABLE))
   {
      outFile << "# histogram is INTEGRATEABLE (bins weights divided by bin width)\n";
   }

   if(norm == Normalize::PDF)
   {
      outFile << "# ... therefore histogram is a probability distribution (PDF)\n";
   }

   outFile << blankLine;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _                                  _ 
    | |/ /  _ \|  _ \| | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___ / |
    | ' /| | | | |_) | |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \| |
    | . \| |_| |  __/|  _  | \__ \ || (_) | (_| | | | (_| | | | | | | |
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_|_|
                                           |___/                       
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

const std::string kdp::KDPHistogram1::FILE_EXTENSION = ".histo1D";
 // An empty set of active names

// Private ctor, for constructing a full object
kdp::KDPHistogram1::KDPHistogram1(std::string const& name,
   weight_t const totalWeight_in,
   weight_t const writeScale_in,
   Normalize const norm_in,
   KDPHistoVec<WeightError>&& disowned):
KDPHistogram_base(name + KDPHistogram1::FILE_EXTENSION,
   totalWeight_in, writeScale_in, norm_in),
hist(std::move(disowned)) {}

// ctor given writeScale
kdp::KDPHistogram1::KDPHistogram1(std::string const& name,
   BinSpecs const& xBinSpecs,
   weight_t const writeScale_in,
   Normalize const norm_in):
KDPHistogram1(name, 0., writeScale_in, norm_in,
   std::move(KDPHistoVec<WeightError>(xBinSpecs))) {}

// ctor given norm
kdp::KDPHistogram1::KDPHistogram1(std::string const& name,
  BinSpecs const& xBinSpecs,
  Normalize const norm_in):
KDPHistogram1(name, xBinSpecs, 1., norm_in) {}


void kdp::KDPHistogram1::Fill(edge_t const x, weight_t const weight,
   weight_t const error)
{
   // x can be infinite, but it cannot be nan
   if(std::isnan(x)) RTExcept("trying to bin nan (x)");
   CheckWeightError(weight, error);
   totalWeight += weight;

   WeightError& bin = hist.FindBin(x);
   bin.weight += weight;
   bin.error2 += error*error;
}

// Divide two histograms to create a new histogram with their ratio
void kdp::KDPHistogram1::WriteRatio(KDPHistogram1 const& that,
   std::string const& newName) const
{
   // Cannot divide histograms unless they have identical BinSpecs
   if(this->hist.binSpecs not_eq that.hist.binSpecs)
      RTExcept("Cannot divide (" + this->fileName + ") by (" +
         that.fileName + ") ... incompatible bin specs.");

   //KDPHistoVec<WeightError> ratio = this->hist;
   auto ratio = this->hist;
   ratio /= that.hist;

   // Create and immediately destroy (and write) the histogram
   // Pass a unit weight (no normalization, so it won't change anything)
   // b/c histo won't auto-write unless it thinks something's there
   KDPHistogram1 ratioHist(newName, 1., 1., Normalize::NO, std::move(ratio));
}

void kdp::KDPHistogram1::Write()
{
   // We already checked we could write to the file,
   // there are no exceptions to generate, and we shouldn't encounter any
   std::ofstream outFile;
   outFile.open(fileName.c_str(), std::ios::out | std::ios::trunc);
   if(outFile)
   {
      BinSpecs const& binSpecs = hist.binSpecs;
      std::vector<edge_t> edges = binSpecs.edges;
      assert(hist.size() == edges.size() - 1);
      
      const bool divideByWeight = bool(norm & Normalize::UNITARY);
      const bool divideByWidth = bool(norm & Normalize::INTEGRATABLE);
      const bool uniform = (binSpecs.binType == BinType::UNIFORM);
      const weight_t normWeight = divideByWeight ? totalWeight : 1.;

      char outBuff[1024]; // should be long enough for one line (ASSUMPTION)
      
      WriteHeader(outFile, "KDPHistogram1");
      
      if(uniform)
      {
         outFile << "# first/last bins are under/over, respectively \n\n";
         // Change the edges we use, so we don't need fancy code
         edges.front() = binSpecs.range.min - binSpecs.binWidth;
         edges.back() = binSpecs.range.max + binSpecs.binWidth;
      }
      else // write under/overflow in the header
      {         
         sprintf(outBuff, "# Under (wgt/err):   %.16e     %.16e \n",
            (writeScale * hist.front().weight) / normWeight,
            (writeScale * std::sqrt(hist.front().error2)) / normWeight);
         outFile << outBuff;

         sprintf(outBuff, "# Over  (wgt/err):   %.16e     %.16e \n\n",
            (writeScale * hist.back().weight) / normWeight,
            (writeScale * std::sqrt(hist.back().error2)) / normWeight);
         outFile << outBuff;
      }

      const int uniformShift = (uniform ? 0 : 1);

      auto const endBin =  hist.end() - uniformShift;
      auto leftEdge = edges.begin() + uniformShift;
      auto rightEdge = leftEdge + 1;      

      outFile << "#     x             pdf(x)        pdfError(x)       cdf(x)         pdf-error       pdf+error\n";

      weight_t cumulative = 0.;      

      for(auto thisBin = hist.begin() + uniformShift; thisBin not_eq endBin; ++thisBin)
      {
         const weight_t thisWeight = thisBin->weight;
         const weight_t thisError = std::sqrt(thisBin->error2);
         const edge_t thisBinDenom = normWeight *
               (divideByWidth ?  (*rightEdge - *leftEdge) : 1.);

         cumulative += thisWeight;

         sprintf(outBuff, "% .6e   % .6e   % .6e   % .6e   % .6e   % .6e\n",
            0.5*(*rightEdge + *leftEdge),
            (writeScale*thisWeight)/thisBinDenom,
            (writeScale*thisError)/thisBinDenom,
            cumulative/totalWeight,
            (writeScale*(thisWeight-thisError))/thisBinDenom,
            (writeScale*(thisWeight+thisError))/thisBinDenom);
         outFile << outBuff;
         
         leftEdge = rightEdge;
         ++rightEdge;
      }

      outFile.close();
      written = true;
   }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     _  ______  ____  _   _ _     _                                  ____  
    | |/ /  _ \|  _ \| | | (_)___| |_ ___   __ _ _ __ __ _ _ __ ___ |___ \ 
    | ' /| | | | |_) | |_| | / __| __/ _ \ / _` | '__/ _` | '_ ` _ \  __) |
    | . \| |_| |  __/|  _  | \__ \ || (_) | (_| | | | (_| | | | | | |/ __/ 
    |_|\_\____/|_|   |_| |_|_|___/\__\___/ \__, |_|  \__,_|_| |_| |_|_____|
                                           |___/                           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 

const std::string kdp::KDPHistogram2::FILE_EXTENSION = ".histo2D";
 // An empty set of active names

// Private ctor, for constructing a full object
kdp::KDPHistogram2::KDPHistogram2(std::string const& name,
   weight_t const totalWeight_in,
   weight_t const writeScale_in,
   Normalize const norm_in,
   KDPHistoVec<KDPHistoVec<WeightError> >&& disowned):
KDPHistogram_base(name + KDPHistogram2::FILE_EXTENSION,
   totalWeight_in, writeScale_in, norm_in),
hist(std::move(disowned)) {}

// Public ctor, through which all other public ctors flow
kdp::KDPHistogram2::KDPHistogram2(std::string const& name,
   BinSpecs const& xBinSpecs, BinSpecs const& yBinSpecs,
   weight_t const writeScale_in,
   Normalize const norm_in):
KDPHistogram2(name, 0., writeScale_in, norm_in,
   std::move(KDPHistoVec<KDPHistoVec<WeightError> >(xBinSpecs, yBinSpecs))) {}

// Minimal ctor (writeScale = 1., norm = PDF)
kdp::KDPHistogram2::KDPHistogram2(std::string const& name,
  BinSpecs const& xBinSpecs, BinSpecs const& yBinSpecs,
  Normalize const norm_in):
KDPHistogram2(name, xBinSpecs, yBinSpecs, 1., norm_in) {}

void kdp::KDPHistogram2::Fill(edge_t const x, edge_t const y,
   weight_t const weight, weight_t const error)
{
   if(std::isnan(x)) RTExcept("trying to bin nan (x)");
   if(std::isnan(y)) RTExcept("trying to bin nan (y)");
   CheckWeightError(weight, error);
   totalWeight += weight;

   WeightError& bin = hist.FindBin(x).FindBin(y);
   bin.weight += weight;
   bin.error2 += error*error;
}

// Divide two histograms to create a new histogram with their ratio
void kdp::KDPHistogram2::WriteRatio(KDPHistogram2 const& that,
   std::string const& newName) const
{
   // Cannot divide histograms unless they have identical BinSpecs
   if((this->hist.binSpecs not_eq that.hist.binSpecs) or
      (this->hist.front().binSpecs not_eq that.hist.front().binSpecs))
      RTExcept("Cannot divide (" + this->fileName + ") by (" +
         that.fileName + ") ... incompatible bin specs.");

   //KDPHistoVec<KDPHistoVec<WeightError> > ratio = this->hist;
   auto ratio = this->hist;
   ratio /= that.hist;

   // Create and immediately destroy (and write) the histogram
   // Pass a unit weight, no normalization or writeScale,
   // but it won't auto-write unless it thinks something's there
   KDPHistogram2 ratioHist(newName, 1., 1., Normalize::NO, std::move(ratio));
}

void kdp::KDPHistogram2::Write()
{
   // We already checked we could write to the file,
   // there are no exceptions to generate, and we shouldn't encounter any
   std::ofstream outFile;
   outFile.open(fileName.c_str(), std::ios::out | std::ios::trunc);
   if(outFile)
   {
      BinSpecs const& binSpecsX = hist.binSpecs;
      BinSpecs const& binSpecsY = hist.front().binSpecs;
      std::vector<edge_t> edgesX = binSpecsX.edges;
      std::vector<edge_t> edgesY = binSpecsY.edges;
      
      const bool divideByWeight = bool(norm & Normalize::UNITARY);
      const bool divideByWidth = bool(norm & Normalize::INTEGRATABLE);

      const bool uniformX = (binSpecsX.binType == BinType::UNIFORM);
      const bool uniformY = (binSpecsY.binType == BinType::UNIFORM);
      
      const weight_t normWeight = divideByWeight ? totalWeight : 1.;

      char outBuff[1024]; // should be long enough for one line (ASSUMPTION)
      
      WriteHeader(outFile, "KDPHistogram2");   

      if(uniformX)
      {
         outFile << "# first/last x-bins are under/over, respectively \n";
         // Change the edges we use, so we don't need fancy code
         edgesX.front() = binSpecsX.range.min - binSpecsX.binWidth;
         edgesX.back() = binSpecsX.range.max + binSpecsX.binWidth;
      }
      else WriteUnderOverX(outFile, normWeight);
      
      if(uniformY)
      {
         outFile << "# first/last y-bins are under/over, respectively \n";
         // Change the edges we use, so we don't need fancy code
         edgesY.front() = binSpecsY.range.min - binSpecsY.binWidth;
         edgesY.back() = binSpecsY.range.max + binSpecsY.binWidth;
      }
      else WriteUnderOverY(outFile, normWeight);

      outFile << "# WARNING: data configured for gnuplot splot or splot with pm3d; \
         the interpolation mesh is probably not being displayed accurately.\n";
      outFile << "# To properly display in gnuplot, you must use \
         \"set pm3d corners2color c1\"\n";
      outFile << "#\n#\n";


      const int uniformShiftX = (uniformX ? 0 : 1);
      const int uniformShiftY = (uniformY ? 0 : 1);

      auto const endHistoX =  hist.end() - uniformShiftX;
      auto leftEdgeX = edgesX.begin() + uniformShiftX;
      auto rightEdgeX = leftEdgeX + 1;

      outFile << "#     x               y             pdf(x)        pdfError(x)\n";

      for(auto thisHistoX = hist.begin() + uniformShiftX;
         thisHistoX not_eq endHistoX; ++thisHistoX)
      {
         //const edge_t xPos = 0.5*(*rightEdgeX + *leftEdgeX);
         // gnuplot assumes you are giving it corners, not central positions
         const edge_t xPos = *leftEdgeX;
         const edge_t xWidth = (*rightEdgeX - *leftEdgeX);
         
         auto const endBinY = thisHistoX->end() - uniformShiftY;
         auto leftEdgeY = edgesY.begin() + uniformShiftY;
         auto rightEdgeY = leftEdgeY + 1;
         
         for(auto thisBinY = thisHistoX->begin() + uniformShiftY;
            thisBinY not_eq endBinY; ++thisBinY)
         {
            const weight_t thisWeight = thisBinY->weight;
            const weight_t thisError = std::sqrt(thisBinY->error2);
            const edge_t thisBinDenom = normWeight *
                  (divideByWidth ?  xWidth*(*rightEdgeY - *leftEdgeY) : 1.);

            sprintf(outBuff, "% .6e   % .6e   % .6e   % .6e\n",
               xPos,
               //0.5*(*rightEdgeY + *leftEdgeY),
               // gnuplot assumes you are giving it corners, not central positions
               *leftEdgeY,
               (writeScale*thisWeight)/thisBinDenom,
               (writeScale*thisError)/thisBinDenom);
            outFile << outBuff;
            
            leftEdgeY = rightEdgeY;
            ++rightEdgeY;
         }
         // gnuplot assumes you are giving it corners, not central positions
         sprintf(outBuff, "% .6e   % .6e   % .6e   % .6e\n",
               xPos,
               //0.5*(*rightEdgeY + *leftEdgeY),
               // gnuplot assumes you are giving it corners, not central positions
               *leftEdgeY, 0., 0.);
            outFile << outBuff;
         
         // gnuplot requires blank line in between x-values
         outFile << "\n";

         leftEdgeX = rightEdgeX;
         ++rightEdgeX;
      }

      // gnuplot assumes you are giving it corners, not central positions
      // Therefore, we must supply the upper boundaries
      for(auto leftEdgeY = edgesY.begin() + uniformShiftY;
         leftEdgeY not_eq edgesY.end() - uniformShiftY; ++leftEdgeY)
      {
         sprintf(outBuff, "% .6e   % .6e   % .6e   % .6e\n",
            *leftEdgeX, *leftEdgeY, 0., 0.);
         outFile << outBuff;
      }

      outFile.close();
      written = true;
   }
}

void kdp::KDPHistogram2::WriteUnderOverX(std::ofstream& outFile,
   weight_t const normWeight)
{
   outFile << "#\n";

   char buffer[1024];

   assert(hist.front().size() == hist.back().size());
   assert(hist.front().size() == hist.front().binSpecs.edges.size() - 1);
   
   auto itUnder = hist.front().begin();
   auto itOver = hist.back().begin();
   auto leftEdge = hist.front().binSpecs.edges.begin();
   auto rightEdge = leftEdge + 1;

   auto const underEnd = hist.front().end();

   while(itUnder not_eq underEnd)
   {
      sprintf(buffer,
         "# x under-over (y-pos | under | over)     %.7e     %.7e     %.7e\n",
         0.5+(*leftEdge + *rightEdge),
         (writeScale * itUnder->weight) / normWeight,
         (writeScale * itOver->weight) / normWeight);
      outFile << buffer;

      ++itUnder;
      ++itOver;
      leftEdge = rightEdge;
      ++rightEdge;
   }
}

void kdp::KDPHistogram2::WriteUnderOverY(std::ofstream& outFile,
   weight_t const normWeight)
{
   outFile << "#\n";

   char buffer[1024];

   // Just assume each vector has the same weight
   assert(hist.size() == hist.binSpecs.edges.size() - 1);
   
   auto leftEdge = hist.binSpecs.edges.begin();
   auto rightEdge = leftEdge + 1;

   for(auto itXhist = hist.begin(); itXhist not_eq hist.end(); ++itXhist)
   {
      sprintf(buffer,
         "# y under-over (x-pos | under | over)     %.7e     %.7e     %.7e\n",
         0.5+(*leftEdge + *rightEdge),
         (writeScale * (itXhist->front().weight)) / normWeight,
         (writeScale * (itXhist->back().weight)) / normWeight);
      outFile << buffer;

      leftEdge = rightEdge;
      ++rightEdge;
   }
}
