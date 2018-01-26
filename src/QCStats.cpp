#include "QCStats.h"
#include "Sequence.h"
#include "QSamFlag.h"
#include "CigarRoller.h"
#include "BaseAsciiMap.h"

#include <cmath>
#include <set>

#define MIN_MAPQ 10
#define SIZE_RESERVED 10000

const int QCStats::depthThreshold[6] = {1, 5, 10, 15, 25, 30};

void QCStats::constructorClear()
{
  for(int i=0;i< 100; i++) {
    coveredGCFreq[i] = 0;
    depthTotalVsGC[i] = 0;
    depthVsGC[i] = 0;
    depthVsGC_norm[i] = 0;
  }
  for(int i=0; i<6; i++) {
    for(int j=0; j<6; j++) {
      matchMatrix[i][j] = 0;
    }
    baseComposition[i] = 0;
  }

}


QCStats::QCStats()
{
  GC=NULL;
  depthVec=NULL;
  constructorClear();
  Init(SIZE_RESERVED);
}

QCStats::QCStats(int n)
{
  if(n<=0) error("Constructor argument has to be >0!\n");
  GC=NULL;
  depthVec=NULL;

  Init(n);
}

QCStats::~QCStats()
{
#if defined(BROKEN)
  // XXX broken - bad memory management, maybe overruning
  // buffers - for now eliminate this and fix the overall
  // program first.
  delete [] matchCountByCycle;
  delete [] misMatchCountByCycle;
  delete [] misMatchRateByCycle;
  delete [] baseCountByQualByCycle;
  for(int i=0; i<size; i++)
  {
    delete [] baseCountByCycle[i];
    delete [] baseCompositionByCycle[i];
  }
  delete [] baseCountByCycle;
  delete [] baseCompositionByCycle;
#endif
}

void QCStats::ReportWarningCount()
{
  if(nWarnings>0) fprintf(stderr, "Total number of warnings: %d \n", (int)nWarnings);
}

void QCStats::Init(int n)
{
  size=0;
  size_reserved = n;
  totalMappedBases = 0;
  matchCountByCycle = new uint64_t[size_reserved];
  misMatchCountByCycle = new uint64_t[size_reserved];
  misMatchRateByCycle = new double[size_reserved];
  baseCountByQualByCycle = new std::map<int, uint64_t>[size_reserved];
  baseQ20CountByCycle.resize(size_reserved);
  baseReportedQ20CountByCycle.resize(size_reserved);

  dbSNP = NULL;

  for(int i=0; i<size_reserved; i++) {
    matchCountByCycle[i] = 0;
    misMatchCountByCycle[i] = 0;
    misMatchRateByCycle[i] = -1.0;
  }
  // GC content
  for(int i=0; i<=100; i++)
  {
    depthTotalVsGC[i] = 0;
    depthVsGC[i] = 0.0;
    depthVsGC_norm[i] = 0.0;
  }

  // Base transition matrix
  for(int i=0; i<6; i++)
    for(int j=0; j<6; j++)
      matchMatrix[i][j] = 0;

  // Base composition along cycles
  baseCountByCycle = new uint64_t *[size_reserved];
  baseCompositionByCycle = new double *[size_reserved];
  for(int i=0; i<size_reserved; i++)
  {
    baseCountByCycle[i] = new uint64_t[6];
    baseCompositionByCycle[i] = new double[6];
  }
  for(int i=0; i<size_reserved; i++)
    for(int j=0; j<6; j++)
    {
      baseCountByCycle[i][j] = 0;
      baseCompositionByCycle[i][j] = 0;
    }

  //For general stats
  nReads=nUnMapped=nUnMapped_Filter=nZeroMapQual=nLT10MapQual=nReadsMapped2TargetRegions=nQ20=nSecondary=nDup=nQCFail=nPaired=nProperPaired=0;
  nBaseCovered  = 0;
  MAX_ISIZE_ALLOWED = INT_MAX; nWarnings = 0;
}

void QCStats::ReAllocateMemory(int size_new)
{
  fprintf(stderr, "Reallocating memroy...\n");
  uint64_t *matchCountByCycle_new = new uint64_t[size_new];
  uint64_t *misMatchCountByCycle_new = new uint64_t[size_new];
  double *misMatchRateByCycle_new = new double[size_new];

  for(int i=0; i<size_new; i++) {
    if(i<size_reserved) {
      matchCountByCycle_new[i] = matchCountByCycle[i];
      misMatchCountByCycle_new[i] = misMatchCountByCycle[i];
      misMatchRateByCycle_new[i] = misMatchRateByCycle[i];
      continue;
    }
    matchCountByCycle_new[i] = 0;
    misMatchCountByCycle_new[i] = 0;
    misMatchRateByCycle_new[i] = -1.0;
  }
  delete [] matchCountByCycle;
  delete [] misMatchCountByCycle;
  delete [] misMatchRateByCycle;
  matchCountByCycle = matchCountByCycle_new;
  misMatchCountByCycle = misMatchCountByCycle_new;
  misMatchRateByCycle = misMatchRateByCycle_new;

  // Base composition along cycles
  uint64_t **baseCountByCycle_new = new uint64_t *[size_new];
  for(int i=0; i<size_new; i++)
    baseCountByCycle_new[i] = new uint64_t[6];
  for(int i=0; i<size_new; i++)
  {
    for(int j=0; j<6; j++)
      if(i<size_reserved)
        baseCountByCycle_new[i][j] = baseCountByCycle[i][j];
      else
        baseCountByCycle_new[i][j] = 0;
  }
  for(int i=0; i<size_reserved; i++)
    delete [] baseCountByCycle[i];
  delete [] baseCountByCycle;
  baseCountByCycle = baseCountByCycle_new;

  size_reserved = size_new;
}

void QCStats::CalcGenomeCoverage(std::vector<bool> &position, uint32_t baseNCount)
{
  nBaseCovered = 0;
  for(uint32_t i=0; i<position.size(); i++)
    if(position[i])
      nBaseCovered++;

  coverage = double(totalMappedBases)/nBaseCovered;
  //coverage = double(totalMappedBases)/(position.size()-baseNCount);
  genomeCoverage = 100*double(nBaseCovered)/(position.size()-baseNCount);
}

void QCStats::CalcMisMatchRateByCycle()
{
  for(int i=0; i<size; i++)
  {
    if((matchCountByCycle[i]+misMatchCountByCycle[i])>0)
      misMatchRateByCycle[i] = double(misMatchCountByCycle[i])/double(matchCountByCycle[i]+misMatchCountByCycle[i]);
  }
}

double QCStats::CalcMisMatchRateByCycle_MEAN()
{
  double sumEPS = 0;
  double eps;
  int cnt=0;
  for(int i=0; i<size; i++)
  {
    if((matchCountByCycle[i]+misMatchCountByCycle[i])>0)
      eps = -10*log10(double(misMatchCountByCycle[i])/double(matchCountByCycle[i]+misMatchCountByCycle[i]));
    else continue;

    if(eps>40) eps = 40;
    sumEPS +=  eps;
    cnt++;
  }
  return(sumEPS/cnt);
}

void QCStats::CalcMisMatchRateByQual()
{
  std::map<int, uint64_t>::iterator p;
  // store all possible qual values
  std::set<int> possibleQual;
  for(p=matchCountByQual.begin(); p!=matchCountByQual.end(); p++)
    possibleQual.insert(p->first);
  for(p=misMatchCountByQual.begin(); p!=misMatchCountByQual.end(); p++)
    possibleQual.insert(p->first);

  // store quality in increasing order to variable "qual"
  qual.clear();
  for (std::set<int>::iterator iter = possibleQual.begin();
       iter != possibleQual.end();
       iter++)
    qual.push_back(*iter);

  for(unsigned int i=0; i<qual.size(); i++)
  {
    qualCount[qual[i]] = matchCountByQual[qual[i]] + misMatchCountByQual[qual[i]];
    if(qualCount[qual[i]]>0) misMatchRateByQual[qual[i]] = double(misMatchCountByQual[qual[i]])/double(qualCount[qual[i]]);
    else misMatchRateByQual[qual[i]] = 0;
  }
}

double QCStats::CalcMisMatchRateByQual_MSE()
{
  if(qual.size()==0) error("No bases are processed!\n");

  double mse = 0;
  uint64_t count = 0;
  for(unsigned int i=0; i<qual.size(); i++)
  {
    if(qual[i]<5 || qualCount[qual[i]]==0) continue;
    double eps = -10*log10(misMatchRateByQual[qual[i]]) > 40 ? 40 : -10*log10(misMatchRateByQual[qual[i]]) ;
    mse += pow((eps - qual[i]), 2) * qualCount[qual[i]];
    count += qualCount[qual[i]];
  }
  return(mse/count);
}

void QCStats::CalcQ20Bases()
{
  uint64_t total = 0;
  nQ20 = 0;
  for(unsigned int i=0; i<qual.size(); i++){
    if(misMatchRateByQual[qual[i]]<=0.01)
    {
      nQ20+=qualCount[qual[i]];
      Q20QualScores[qual[i]]++;
    }
    total += qualCount[qual[i]];
  }
  pQ20 = 100*double(nQ20)/total;
}

void QCStats::CalcQ20BasesByCycle()
{
  for(int cycle=0; cycle<size; cycle++)
  {
    std::vector<int> qualInACycle;
    std::map<int, uint64_t>::iterator p;

    for(p=baseCountByQualByCycle[cycle].begin(); p!=baseCountByQualByCycle[cycle].end(); p++)
      qualInACycle.push_back(p->first);

    for(unsigned int i=0; i<qualInACycle.size(); i++)
    {
      if(Q20QualScores[qualInACycle[i]]>0)
        baseQ20CountByCycle[cycle] += baseCountByQualByCycle[cycle][qualInACycle[i]];
    }
  }
}

void QCStats::CalcReportedQ20BasesByCycle()
{
  for(int cycle=0; cycle<size; cycle++) {
    if(baseReportedQ20CountByCycle[cycle]>0) error("Reported Q20 bases have been calculated!\n");
  }

  for(int cycle=0; cycle<size; cycle++)
  {
    std::vector<int> qualInACycle;
    std::map<int, uint64_t>::iterator p;

    for(p=baseCountByQualByCycle[cycle].begin(); p!=baseCountByQualByCycle[cycle].end(); p++)
      qualInACycle.push_back(p->first);

    for(unsigned int i=0; i<qualInACycle.size(); i++)
    {
      if(qualInACycle[i]>=20)
        baseReportedQ20CountByCycle[cycle] += baseCountByQualByCycle[cycle][qualInACycle[i]];
    }
  }
}

void QCStats::CalcBaseComposition()
{
  uint64_t total = 0;
  uint64_t total_byCycle = 0;
  for(int i=0; i<size; i++)
  {
    total_byCycle = 0;
    for(int j=0; j<6; j++)
    {
      baseComposition[j] += baseCountByCycle[i][j];
      total+=baseCountByCycle[i][j];
      total_byCycle += baseCountByCycle[i][j];
    }
    for(int j=0; j<6; j++)  {
      if(total_byCycle>0)	baseCompositionByCycle[i][j] = baseCountByCycle[i][j]/total_byCycle;
      else	baseCompositionByCycle[i][j] = 0;
    }
  }
  for(int j=0; j<6; j++)
    if(total>0)
      baseComposition[j] = 100*double(baseComposition[j])/total;
    else baseComposition[j] = 0;
}

void QCStats::CalcDepthGC(GCContent &gc, std::vector<bool> &genomePosCovered)
{
  for(int i=0; i<=100; i++) coveredGCFreq[i] = 0;
  for(uint32_t i=0; i<genomePosCovered.size(); i++)
    if(genomePosCovered[i])
      coveredGCFreq[gc.gcCount[i]]++;

  for(int i=0; i<=100; i++)
    if(coveredGCFreq[i]>0){
      depthVsGC[i] = double(depthTotalVsGC[i])/coveredGCFreq[i];
      depthVsGC_norm[i] = depthVsGC[i]/coverage;
    }
}

void QCStats::CalcDepthDist()
{
  // reset depth distribution counter
  int nThreshold = sizeof(depthThreshold)/sizeof(depthThreshold[0]);
  for(int i = 0; i < nThreshold; ++i) {
    depthDistribution[i] = 0;
  }
  
  const std::vector<uint32_t>& freq = depthVec->getFreqDist();
  for (unsigned int i = 1; i < freq.size(); i++){
    depthDist[i] = freq[i];
    // fprintf(stderr, "depth %u count = %u\n", i, freq[i]);

    for(int j = 0; j < nThreshold; ++j) {
      if (i >= (unsigned int)depthThreshold[j]) {
        depthDistribution[j] += freq[i];
      }
    }
  }

  // for(int i = 0; i < nThreshold; ++i) {
  //   fprintf(stderr, "depth > %d sites: %d\n", depthThreshold[i], depthDistribution[i]);
  // }
}

//Should be called after CalcDepthGC and CalcDepthDist to get relative info
void QCStats::CalcGCBias(int start, int end)
{
  uint64_t totalCoveredGC = 0;
  double var = 0.0;
  for(int i=start; i<=end; i++)
  {
    var += pow((depthVsGC_norm[i]-1),2)* coveredGCFreq[i];
    totalCoveredGC += coveredGCFreq[i];
  }
  gcBiasStat = var/totalCoveredGC;
}

void QCStats::CalcInsertSize_mean()
{
  uint64_t totalInsertSize = 0;
  uint64_t totalInsertCount = 0;
  std::map<int32_t, uint64_t>::iterator p;
  for(p=insertSize.begin(); p!=insertSize.end();p++)
  {
    totalInsertSize += (p->first * p->second);
    totalInsertCount += (p->second);
  }
  insertSize_mean = double(totalInsertSize)/totalInsertCount;
}

void QCStats::CalcInsertSize_var()
{
  uint64_t totalSSE = 0;
  uint64_t totalInsertCount = 0;
  std::map<int32_t, uint64_t>::iterator p;
  for(p=insertSize.begin(); p!=insertSize.end();p++)
  {
    totalSSE += pow((p->first-insertSize_mean),2) * p->second;
    totalInsertCount += p->second;
  }
  insertSize_var = double(totalSSE)/totalInsertCount;
}

void QCStats::CalcInsertSize_mode()
{
  uint64_t maxCnt = 0;
  insertSize_mode = 0;
  std::map<int32_t, uint64_t>::iterator p;
  for(p=insertSize.begin(); p!=insertSize.end();p++)
    if(p->second>maxCnt)
    {
      maxCnt = p->second;
      insertSize_mode = p->first;
    }
}

void QCStats::CalcInsertSize_medium()
{
  uint64_t totalCnt = 0;
  uint64_t partialCnt = 0;
  insertSize_medium = 0;

  std::map<int32_t, uint64_t>::iterator p;
  for(p=insertSize.begin(); p!=insertSize.end();p++)
    totalCnt+=p->second;

  for(p=insertSize.begin(); p!=insertSize.end();p++)
  {
    partialCnt += p->second;
    if(partialCnt>=totalCnt/2)  {
      insertSize_medium = p->first;
      return;
    }
  }
}

void QCStats::PrintSamRecord(SamRecord &sam)
{
  fprintf(stderr, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
          (const char *) sam.getReadName(),
          (int)sam.getFlag(),
          (const char *) sam.getReferenceName(),
          sam.get1BasedPosition(),
          sam.getMapQuality(),
          (const char *) sam.getCigar(),
          (const char *) sam.getMateReferenceName(),
          sam.get1BasedMatePosition(),
          sam.getInsertSize(),
          (const char *) sam.getSequence(),
          (const char *) sam.getQuality());

}

void QCStats::UpdateStats(SamRecord & sam, QSamFlag &filter, double minMapQuality, std::map<int, int> &lanes, std::vector<std::string>& readGroups)
{
  if(lanes.size()>0) {
    StringArray tokens;
    tokens.Clear();
    tokens.AddTokens(sam.getReadName(), ":");
    if(lanes[tokens[1].AsInteger()]==0) return;
  };

  if (readGroups.size() > 0) {
    const String* p_rg = sam.getStringTag("RG");
    if (p_rg == NULL) return;
    const String& rg = *p_rg;
    bool match = true;
    for (unsigned int i = 0; i < readGroups.size(); i++) {
      for (unsigned int j = 0; j < readGroups[i].size(); j++ ){
        if (j == (unsigned int) rg.Length()) {
          match = false;
          break;
        }
        if (readGroups[i][j] != rg[j]) {
          match = false;
          break;
        }
      }
      if (match) break;
    }
    if (!match) return;
  };

  QSamFlag flag;
  flag.GetFlagFields(sam.getFlag());


  // Filters
  if(flag.isPaired==true)
  {
    if(filter.isRead1==true && flag.isRead1==true) return;
    if(filter.isRead2==true && flag.isRead2==true) return;
  }
  if(filter.isPaired==true && flag.isPaired==true) return;
  if(filter.isUnPaired==true && flag.isPaired==false) return;

  nReads++;

  String refLabel = sam.getReferenceName();
  //  genomeIndex_t mapPos = referencegenome->getGenomePosition(refLabel.c_str(), sam.data.header->position+1);
  genomeIndex_t mapPos = referencegenome->getGenomePosition(refLabel.c_str(), sam.get1BasedPosition());

  if(sam.getReadLength()>size){
    if(sam.getReadLength()>size_reserved)
      ReAllocateMemory(sam.getReadLength()*2);
    size = sam.getReadLength();
  }

  if(flag.isPaired) {
    nPaired++;
    if(flag.isProperPaired)
      nProperPaired++;
  }
  if(flag.isUnMapped) { nUnMapped++; nUnMapped_Filter++; }
  else {
    if(sam.getMapQuality()==0) nZeroMapQual++;
    if(sam.getMapQuality()<10) nLT10MapQual++;
  }
  if(flag.isSecondary) { nSecondary++;  }
  if(flag.isDup) { nDup++;  }
  if(flag.isQCFail) { nQCFail++;  }

  if(flag.isUnMapped) { return; }
  if(filter.isSecondary && flag.isSecondary) { return;  }
  if(filter.isDup && flag.isDup) { return;  }
  if(filter.isQCFail && flag.isQCFail) { return;  }

  if(sam.getMapQuality()<minMapQuality) { nUnMapped_Filter++; return; }

  // Update insert size for paired reads
  if(insertSize.size()==10000) {
    //CalcInsertSize_medium();
    std::map<int32_t, uint64_t>::iterator p;
    p = insertSize.end();
    MAX_ISIZE_ALLOWED = (--p)->first > 10000? (--p)->first : 10000;
  }

  if(flag.isPaired)
  {
    int iSize = sam.getInsertSize();
    if(iSize > MAX_ISIZE_ALLOWED) iSize =  MAX_ISIZE_ALLOWED;
    if(sam.getMapQuality()>MIN_MAPQ && iSize>0)
      insertSize[iSize]++;
  }


  // Update cycle number since different lenght of read may exist in the same bam file
  // Softclips are excluded from the read length
  //  cycles[sam.data.sequence.Length()-cigarcol.nClips_begin-cigarcol.nClips_end]++;
  cycles[sam.getReadLength()]++;

  if(mapPos==INVALID_GENOME_INDEX)
  {
    if(++nWarnings<=5)
    {
      fprintf(stderr, "WARNING: INVALID_GENOME_INDEX (%u) and record skipped... Reference in BAM is different from the ref used here!\n", mapPos);
      PrintSamRecord(sam);
      if(nWarnings==5) fprintf(stderr, "Too many warnings and rest of them are not reported!\n");
    }
    return;
  }

  // Get the cigar.
  Cigar* cigar = sam.getCigarInfo();

  String aligTypes="";
  for (int i = 0; i < cigar->size(); i++)
    for (unsigned int j = 0; j < cigar->getOperator(i).count; j++)
      aligTypes += cigar->getOperator(i).getChar();

  int offset = 0; //to adjust position in ref genome due to insertion 'I' and 'S'
  int offset2 = 0;  //to adjust position in a read due to deletion 'D'
  int offset3 = 0; //to adjust positions due to hard clip 'H' at the start
  int offset4 = 0; //to adjust positions due to hard clip 'H' at the end
  int nCycles = sam.getReadLength();

  char refBase, readBase;
  uint32_t refpos;
  int seqpos;
  int cycleIdx;
  int qual;
  bool coverRegion = false;

  bool isNewRead = true;
  bool hardclip_head = true;

  for (int i = 0; i < aligTypes.Length(); i++) {
    if (aligTypes[i] == 'S') {
      offset--;
      hardclip_head = false;
      continue;
    }
    else if (aligTypes[i] == 'I') {
      offset--;
      hardclip_head = false;

      continue;
    }
    else if (aligTypes[i] == 'D' || aligTypes[i] == 'N') {
      offset2--;
      hardclip_head = false;
      continue;
    }
	else if(aligTypes[i] == 'H') {
      if(hardclip_head) offset3--;
      else offset4--;
      continue;
	}

    refpos = mapPos + i + offset3 + offset;

    if((*regionIndicator).size()>0)
    {
      if((*regionIndicator)[refpos]==false) continue;
      coverRegion = true;
    }

    // Update depth vector and GC content vector
    if (isNewRead){
      isNewRead = false;
      depthVec->beginNewRead(refpos);
    }
    depthVec->addBase(refpos);
    // if(depthVec!=NULL)
    //   if((*depthVec).depth[refpos]<255)(*depthVec).depth[refpos]++;
    if(GC!=NULL)
      depthTotalVsGC[GC->gcCount[refpos]]++;

    (*genomePosCovered)[refpos]=true;
    if(refBase!='N') totalMappedBases++;
    // Excluding dbSNPs for mismatch rate calculation
    if((*dbSNP).size()>0 && (*dbSNP)[refpos]==true) continue;

    refBase = toupper((*referencegenome)[refpos]);
    if(flag.isReverse==true) refBase = BaseAsciiMap::base2complement[(uint32_t) refBase];

    seqpos = i + offset2 + offset3;
    readBase = toupper(sam.getSequence()[seqpos]);

    cycleIdx = seqpos - offset3;
    //printf("CIGAR=%s:Len=%d\nrefpos=%d offset=%d offset2=%d offset3=%d offset4=%d seqpos=%d cycleIdx=%d isReverse=%d\n", aligTypes.c_str(), aligTypes.Length(), refpos, offset, offset2, offset3, offset4, seqpos, cycleIdx, flag.isReverse);

    if(flag.isReverse==true){
      cycleIdx = nCycles - offset3 - offset4 - cycleIdx - 1;
      readBase = BaseAsciiMap::base2complement[(uint32_t) readBase];
    }

    // Maped based
    //if(refBase=='A' || refBase=='C' || refBase=='G' || refBase=='T') totalMappedBases++;
//    if(refBase!='N') totalMappedBases++;
    //
    // Base compostion
    //
    // NB: when bases match the reference, they are returned as '=',
    // which base2int maps to 5.  But the libsrc Sam code returns '0'
    // instead, which maps to 0 (because BaseAsciiMap::base2int
    // maps color space values as well as base space values).
    //
    if(readBase!='=' && readBase != '0') baseCountByCycle[cycleIdx][BaseAsciiMap::base2int[(uint32_t) readBase]]++;
    // Base transition matrix (not yet implemented the plotting)
    // matchMatrix[Sequence::nt2idx(refBase)][Sequence::nt2idx(readBase)]++;

    qual = sam.getQuality()[seqpos]-33;

    //update base count for each quality on each cycle
    baseCountByQualByCycle[cycleIdx][qual]++;

    // update match/mismatch counts
    if(refBase==readBase || readBase=='=')
    {
      matchCountByCycle[cycleIdx]++;
      matchCountByQual[qual]++;
    }
    else if(refBase!='N' && readBase!='N'){
      misMatchCountByCycle[cycleIdx]++;
      misMatchCountByQual[qual]++;
    }
  }
  if((*regionIndicator).size()>0 && coverRegion==true) nReadsMapped2TargetRegions++;
}
