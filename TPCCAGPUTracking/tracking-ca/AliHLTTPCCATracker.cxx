// @(#) $Id$
//***************************************************************************
// This file is property of and copyright by the ALICE HLT Project          * 
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  Ivan Kisel <kisel@kip.uni-heidelberg.de>                *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//***************************************************************************

#include "AliHLTTPCCATracker.h"
#include "AliHLTTPCCAOutTrack.h"
#include "AliHLTTPCCAGrid.h"
#include "AliHLTTPCCARow.h"
#include "AliHLTTPCCATrack.h"
#include "AliHLTTPCCATracklet.h"
#include "AliHLTTPCCAMath.h"
#include "AliHLTTPCCAHit.h"

#include "TStopwatch.h"
#include "AliHLTTPCCAHitArea.h"
#include "AliHLTTPCCANeighboursFinder.h"
#include "AliHLTTPCCANeighboursCleaner.h"
#include "AliHLTTPCCAStartHitsFinder.h"
#include "AliHLTTPCCATrackletConstructor.h"
#include "AliHLTTPCCATrackletSelector.h"
#include "AliHLTTPCCAProcess.h"
#include "AliHLTTPCCAUsedHitsInitialiser.h"

#include "AliHLTTPCCATrackParam.h"

#if !defined(HLTCA_GPUCODE)
#include <iostream>
#endif

//#define DRAW

#ifdef DRAW
  #include "AliHLTTPCCADisplay.h"
  #include "TApplication.h"
#endif //DRAW

ClassImp(AliHLTTPCCATracker)

#if !defined(HLTCA_GPUCODE)  

AliHLTTPCCATracker::AliHLTTPCCATracker()
  :
  fParam(),
  fNHitsTotal(0),  
  fCommonMemory(0), 
  fCommonMemorySize(0),
  fInputEvent(0),     
  fInputEventSize(0), 
  fRowData(0),     
  fRowDataSize(0), 
  fHitInputIDs(0), 
  fHitWeights(0),  
  fNTracklets(0),
  fTrackletStartHits(0),
  fTracklets(0), 
  fNTracks(0),   
  fTracks(0), 
  fNTrackHits(0),
  fTrackHits(0),
  fNOutTracks(0),
  fOutTracks(0), 
  fNOutTrackHits(0),
  fOutTrackHits(0)
{
  // constructor
}

AliHLTTPCCATracker::AliHLTTPCCATracker( const AliHLTTPCCATracker& )
  :
  fParam(),
  fNHitsTotal(0),
  fCommonMemory(0), 
  fCommonMemorySize(0),
  fInputEvent(0),     
  fInputEventSize(0), 
  fRowData(0),     
  fRowDataSize(0), 
  fHitInputIDs(0), 
  fHitWeights(0),  
  fNTracklets(0),
  fTrackletStartHits(0),
  fTracklets(0), 
  fNTracks(0),   
  fTracks(0), 
  fNTrackHits(0),
  fTrackHits(0),
  fNOutTracks(0),
  fOutTracks(0), 
  fNOutTrackHits(0),
  fOutTrackHits(0)
{
  // dummy
}

AliHLTTPCCATracker &AliHLTTPCCATracker::operator=( const AliHLTTPCCATracker& )
{
  // dummy
  fCommonMemory = 0;
  return *this;
}

GPUd() AliHLTTPCCATracker::~AliHLTTPCCATracker()
{
  // destructor
  if( fCommonMemory ) delete[] fCommonMemory;
}
#endif



// ----------------------------------------------------------------------------------
GPUd() void AliHLTTPCCATracker::Initialize( AliHLTTPCCAParam &param )
{
  // initialisation
  fParam = param;
  fParam.Update();
  for( Int_t irow=0; irow<fParam.NRows(); irow++ ){
    fRows[irow].X() = fParam.RowX(irow);
    fRows[irow].MaxY() = CAMath::Tan( fParam.DAlpha()/2.)*fRows[irow].X();
  }
  StartEvent();
}

GPUd() void AliHLTTPCCATracker::StartEvent()
{
  // start new event and fresh the memory  

  fNHitsTotal = 0;
  Int_t oldSize = fCommonMemorySize;
  SetPointers();
  if( fCommonMemory && fCommonMemorySize > oldSize ){
    delete[] fCommonMemory;
    fCommonMemory = 0;
  }
  if( !fCommonMemory ){   
    fCommonMemory = reinterpret_cast<Char_t*> ( new uint4 [ fCommonMemorySize/sizeof(uint4) + 100] );
  }
  SetPointers();  
  *fNTracklets = 0;
  *fNTracks = 0 ;
  *fNTrackHits = 0;
  *fNOutTracks = 0;
  *fNOutTrackHits = 0;
}

GPUhd() void  AliHLTTPCCATracker::SetPointers()
{
  // set all pointers to the event memory

  Int_t gridSizeTotal = 2*(2*fNHitsTotal + 10*Param().NRows());

  ULong_t mem = (ULong_t) fCommonMemory;  
  UInt_t sI = sizeof(Int_t);
  UInt_t sF = sizeof(Float_t);
  UInt_t sS = sizeof(Short_t);
  UInt_t s4 = sizeof(uint4);

  fInputEvent = (Char_t*) mem;  
  fInputEventSize = (1+fParam.NRows()*2 + 1)*sI + (fNHitsTotal*2)*sF;
  mem+= fInputEventSize;

  mem = ( mem/s4 + 1 )*s4;
  fRowData = (uint4*) mem;
  fRowDataSize = ( 2*fNHitsTotal*sS +  //  yz
		   gridSizeTotal*sS + // grid
		   2*fNHitsTotal*sS +  // link up,link down
		   fParam.NRows()*s4   // row alignment
		   );
  mem += fRowDataSize;

  mem = ( mem/sI + 1 )*sI;

  fHitInputIDs = (Int_t*) mem;
  mem+= fNHitsTotal*sI;

  fHitWeights = (Int_t*) mem;
  mem+=  fNHitsTotal*sI;

  fNTracklets = (Int_t*) mem;
  mem+= sI;

  fTrackletStartHits = (Int_t*) mem;
  mem+= fNHitsTotal*sI;

  mem = ( mem/sizeof(AliHLTTPCCATracklet) + 1 )*sizeof(AliHLTTPCCATracklet);

  fTracklets = (AliHLTTPCCATracklet *) mem;
  mem+= fNHitsTotal*sizeof(AliHLTTPCCATracklet);
  
  mem = ( mem/sI + 1 )*sI;

  fNTracks = (Int_t*) mem;
  mem+= sI;

  mem = ( mem/sizeof(AliHLTTPCCATrack) + 1 )*sizeof(AliHLTTPCCATrack);

  fTracks = (AliHLTTPCCATrack*) mem;
  mem+= fNHitsTotal*sizeof(AliHLTTPCCATrack);
  
  mem = ( mem/sI + 1 )*sI;
  fNTrackHits = (Int_t*) mem;
  mem+= sI;

  fTrackHits = (Int_t*) mem;
  mem+= fNHitsTotal*sI;

  fNOutTracks = (Int_t*) mem;
  mem+= sI;

  mem = ( mem/sizeof(AliHLTTPCCAOutTrack) + 1 )*sizeof(AliHLTTPCCAOutTrack);
  
  fOutTracks = (AliHLTTPCCAOutTrack*) mem;
  mem+= fNHitsTotal*sizeof(AliHLTTPCCAOutTrack);

  mem = ( mem/sI + 1 )*sI;
  fNOutTrackHits = (Int_t*) mem;
  mem+= sI;
  fOutTrackHits = (Int_t*) mem;
  mem+= fNHitsTotal*sI;

  fCommonMemorySize = mem - (ULong_t) fCommonMemory;
}


GPUd() void AliHLTTPCCATracker::ReadEvent( Int_t *RowFirstHit, Int_t *RowNHits, Float_t *Y, Float_t *Z, Int_t NHits )
{
  //* Read event

  fNHitsTotal = NHits;

  {
    Int_t oldSize = fCommonMemorySize;
    SetPointers();
    if( fCommonMemory && fCommonMemorySize > oldSize ){
      delete[] fCommonMemory;
      fCommonMemory = 0;
    }
    if( !fCommonMemory ){   
      fCommonMemory = reinterpret_cast<Char_t*> ( new uint4 [ fCommonMemorySize/sizeof(uint4) + 100] );
    }
    SetPointers();  
    *fNTracklets = 0;
    *fNTracks = 0 ;
    *fNOutTracks = 0;
    *fNOutTrackHits = 0;
  }

  reinterpret_cast<Int_t*>( fInputEvent )[0] = fParam.NRows();
  reinterpret_cast<Int_t*>( fInputEvent )[1+fParam.NRows()*2] = NHits;
  Int_t *rowHeaders = reinterpret_cast<Int_t*>( fInputEvent ) +1;
  Float_t *hitsYZ = reinterpret_cast<Float_t*>( fInputEvent ) + 1+fParam.NRows()*2+1;
  for( Int_t iRow=0; iRow<fParam.NRows(); iRow++ ){
    rowHeaders[iRow*2  ] = RowFirstHit[iRow];
    rowHeaders[iRow*2+1] = RowNHits[iRow];
  }
  for( Int_t iHit=0; iHit<NHits; iHit++ ){
    hitsYZ[iHit*2  ] = Y[iHit];
    hitsYZ[iHit*2+1] = Z[iHit];
  }
  SetupRowData();
}


GPUd() void AliHLTTPCCATracker::SetupRowData()
{
  //* Convert input hits, create grids, etc.

  fNHitsTotal = reinterpret_cast<Int_t*>( fInputEvent )[1+fParam.NRows()*2];
  Int_t *rowHeaders = reinterpret_cast<Int_t*>( fInputEvent ) +1;
  Float_t *hitsYZ = reinterpret_cast<Float_t*>( fInputEvent ) + 1+fParam.NRows()*2+1;
  for( Int_t iRow=0; iRow<fParam.NRows(); iRow++ ){
    AliHLTTPCCARow &row = fRows[iRow];
    row.FirstHit() = rowHeaders[iRow*2];
    row.NHits() = rowHeaders[iRow*2+1];
    Float_t yMin=1.e3, yMax=-1.e3, zMin=1.e3, zMax=-1.e3;
    Int_t nGrid =  row.NHits();
    for( Int_t i=0; i<row.NHits(); i++ ){
      Int_t j = row.FirstHit()+i;
      Float_t y = hitsYZ[j*2];
      Float_t z = hitsYZ[j*2+1];
      if( yMax < y ) yMax = y;
      if( yMin > y ) yMin = y;
      if( zMax < z ) zMax = z;
      if( zMin > z ) zMin = z;
    }
    if( nGrid <= 0 ){
      yMin = yMax = zMin = zMax = 0;
      nGrid = 1;
    }

    row.Grid().Create( yMin, yMax, zMin, zMax, nGrid );
    
    Float_t sy = ( CAMath::Abs( row.Grid().StepYInv() ) >1.e-4 ) ?1./row.Grid().StepYInv() :1;
    Float_t sz = ( CAMath::Abs( row.Grid().StepZInv() ) >1.e-4 ) ?1./row.Grid().StepZInv() :1;
    
    //cout<<"grid n = "<<row.Grid().N()<<" "<<sy<<" "<<sz<<" "<<yMin<<" "<<yMax<<" "<<zMin<<" "<<zMax<<endl;
    
    Bool_t recreate=0;
    if( sy < 2. ) { recreate = 1; sy = 2; }
    if( sz < 2. ) { recreate = 1; sz = 2; }
    if( recreate ) row.Grid().Create( yMin, yMax, zMin, zMax, sy, sz );
    
    //cout<<"grid row "<<iRow<<", n = "<<row.Grid().N()<<endl;
  }    

  AliHLTTPCCAHit ffHits[fNHitsTotal];

  Int_t rowDataOffset = 0;

  for( Int_t iRow=0; iRow<fParam.NRows(); iRow++ ){
    AliHLTTPCCARow &row = fRows[iRow];
    AliHLTTPCCAGrid &grid = row.Grid();

    Int_t c[grid.N()+3+10];
    Int_t bins[row.NHits()];
    Int_t filled[ row.Grid().N() +3+10 ];

    for( UInt_t bin=0; bin<row.Grid().N()+3; bin++ ) filled[bin] = 0;  

    for( Int_t i=0; i<row.NHits(); i++ ){
      Int_t j = row.FirstHit()+i;
      Int_t bin = row.Grid().GetBin( hitsYZ[2*j], hitsYZ[2*j+1] );
      bins[i] = bin;
      filled[bin]++;
    }

    {
      Int_t n=0;
      for( UInt_t bin=0; bin<row.Grid().N()+3; bin++ ){
	c[bin] = n;
	n+=filled[bin];
      }
    }
    for( Int_t i=0; i<row.NHits(); i++ ){ 
      Int_t bin = bins[i];
      Int_t ind = c[bin] + filled[bin]-1;
      AliHLTTPCCAHit &h = ffHits[row.FirstHit()+ind];
      fHitInputIDs[row.FirstHit()+ind] = row.FirstHit()+i;
      h.Y() = hitsYZ[2*(row.FirstHit()+i)];
      h.Z() = hitsYZ[2*(row.FirstHit()+i)+1];
      filled[bin]--;
    }

    {
      Float_t y0 = row.Grid().YMin();
      Float_t stepY = (row.Grid().YMax() - y0)*(1./65535.);
      Float_t z0 = row.Grid().ZMin();
      Float_t stepZ = (row.Grid().ZMax() - z0)*(1./65535.);
      Float_t stepYi = 1./stepY;
      Float_t stepZi = 1./stepZ;
      
      row.Hy0() = y0;
      row.Hz0() = z0;
      row.HstepY() = stepY;
      row.HstepZ() = stepZ;
      row.HstepYi() = stepYi;
      row.HstepZi() = stepZi;

      row.FullOffset() = rowDataOffset;
      ushort2 *p= (ushort2*)( fRowData + row.FullOffset() );
      for( Int_t ih=0; ih<row.NHits(); ih++ ){
  	Int_t ihTot = row.FirstHit()+ih;
	AliHLTTPCCAHit &hh = ffHits[ihTot];
	Float_t xx = ((hh.Y() - y0)*stepYi);
	Float_t yy = ((hh.Z() - z0)*stepZi);
	if( xx<0 || yy<0 || xx>=65536 || yy>= 65536 ){
	  std::cout<<"!!!! hit packing error!!! "<<xx<<" "<<yy<<" "<<std::endl;
	}
	p[ih].x = (UShort_t) xx;
	p[ih].y = (UShort_t) yy;
      }
      Int_t size = row.NHits()*sizeof(ushort2);

      row.FullGridOffset() = row.NHits()*2;      
      UShort_t *p1 = ((UShort_t *)p) + row.FullGridOffset();

      Int_t n = grid.N();
      for( Int_t i=0; i<n; i++ ){
	p1[i] = c[i];
      }     
      UShort_t a = c[n];
      Int_t nn = n+grid.Ny()+3;
      for( Int_t i=n; i<nn; i++ ) p1[i] = a;

      size+= (nn)*sizeof(UShort_t);
      row.FullLinkOffset() = row.NHits()*2 + nn;
      size+= row.NHits()*2*sizeof(Short_t);
      if( size%16 ) size = size/sizeof(uint4)+1;
      else size = size/sizeof(uint4);
      row.FullSize()=size;
      //cout<<iRow<<", "<<row.fNHits<<"= "<<size*16<<"b: "<<row.fFullOffset<<" "<<row.fFullSize<<" "<<row.fFullGridOffset<<" "<<row.fFullLinkOffset<<std::endl;

      rowDataOffset+=size;
    }
  } 
}


GPUh() void AliHLTTPCCATracker::Reconstruct()
{
  //* reconstruction of event
   
#ifdef DRAW
  if( !gApplication ){
    TApplication *myapp = new TApplication("myapp",0,0);
  }
  //AliHLTTPCCADisplay::Instance().Init();
  
  AliHLTTPCCADisplay::Instance().SetCurrentSlice( this );
  AliHLTTPCCADisplay::Instance().SetSliceView();
  AliHLTTPCCADisplay::Instance().DrawSlice( this );  
  //for( Int_t iRow=0; iRow<fParam.NRows(); iRow++ )
  //for (Int_t i = 0; i<fRows[iRow].NHits(); i++) 
  //AliHLTTPCCADisplay::Instance().DrawHit( iRow, i );
  //AliHLTTPCCADisplay::Instance().Ask();
#endif

  fTimers[0] = 0; // find neighbours
  fTimers[1] = 0; // construct tracklets
  fTimers[2] = 0; // fit tracklets
  fTimers[3] = 0; // prolongation of tracklets
  fTimers[4] = 0; // selection
  fTimers[5] = 0; // write output
  fTimers[6] = 0;
  fTimers[7] = 0;

  if( fNHitsTotal < 1 ) return;
  TStopwatch timer0;

  *fNTracks = 0;
  *fNTracklets = 0;

#if !defined(HLTCA_GPUCODE)  

  AliHLTTPCCAProcess<AliHLTTPCCANeighboursFinder>( Param().NRows(), 1, *this );
  AliHLTTPCCAProcess<AliHLTTPCCANeighboursCleaner>( Param().NRows()-2, 1, *this );
  AliHLTTPCCAProcess<AliHLTTPCCAStartHitsFinder>( Param().NRows()-4, 1, *this );

  Int_t nStartHits = *fNTracklets;
  
  Int_t nThreads = 128;
  Int_t nBlocks = fNHitsTotal/nThreads + 1;
  if( nBlocks<12 ){
    nBlocks = 12; 
    nThreads = fNHitsTotal/12+1;
    if( nThreads%32 ) nThreads = (nThreads/32+1)*32;
  }
      
  nThreads = fNHitsTotal;
  nBlocks = 1;

  AliHLTTPCCAProcess<AliHLTTPCCAUsedHitsInitialiser>(nBlocks, nThreads,*this);

  Int_t nMemThreads = AliHLTTPCCATrackletConstructor::NMemThreads();
  nThreads = 256;//96;
  nBlocks = nStartHits/nThreads + 1;
  if( nBlocks<30 ){
    nBlocks = 30;
    nThreads = (nStartHits)/30+1;
    if( nThreads%32 ) nThreads = (nThreads/32+1)*32;
  }

  nThreads = nStartHits;
  nBlocks = 1;

  AliHLTTPCCAProcess1<AliHLTTPCCATrackletConstructor>(nBlocks, nMemThreads+nThreads,*this);

  { 
    nThreads = 128;
    nBlocks = nStartHits/nThreads + 1;
    if( nBlocks<12 ){
      nBlocks = 12;  
      nThreads = nStartHits/12+1;
      nThreads = (nThreads/32+1)*32;
    }
    
    *fNTrackHits = 0;

    nThreads = nStartHits;
    nBlocks = 1;


    AliHLTTPCCAProcess<AliHLTTPCCATrackletSelector>(nBlocks, nThreads,*this);

    //cout<<"Slice "<<Param().ISlice()<<": N start hits/tracklets/tracks = "<<nStartHits<<" "<<nStartHits<<" "<<*fNTracks<<std::endl;
    WriteOutput();      
  }

#endif

  timer0.Stop();  
  fTimers[0] = timer0.CpuTime()/100.;

 }




GPUh() void AliHLTTPCCATracker::WriteOutput()
{
  // write output

  TStopwatch timer;
 
  *fNOutTrackHits = 0;
  *fNOutTracks = 0;

  //cout<<"output: nTracks = "<<*fNTracks<<", nHitsTotal="<<fNHitsTotal<<std::endl;  

  for( Int_t iTr=0; iTr<*fNTracks; iTr++){

    AliHLTTPCCATrack &iTrack = fTracks[iTr];

    //cout<<"iTr = "<<iTr<<", nHits="<<iTrack.NHits()<<std::endl;
 
    if( !iTrack.Alive() ) continue;
    if( iTrack.NHits()<3 ) continue;          
    AliHLTTPCCAOutTrack &out = fOutTracks[*fNOutTracks];
    out.FirstHitRef() = *fNOutTrackHits;
    out.NHits() = 0;
    out.OrigTrackID() = iTr;
    {
      out.StartPoint() = iTrack.Param();
      out.EndPoint() = iTrack.Param();
    }

    Int_t iID = iTrack.FirstHitID();
    Int_t nOutTrackHitsOld = *fNOutTrackHits;
 
    for( Int_t ith=0; ith<iTrack.NHits(); ith++ ){
      Int_t ic = (fTrackHits[iID+ith]);
      AliHLTTPCCARow &row = ID2Row(ic);
      Int_t ih = ID2IHit(ic);
      fOutTrackHits[*fNOutTrackHits] = fHitInputIDs[row.FirstHit()+ih];      
      (*fNOutTrackHits)++;
      //cout<<"write i,row,hit,id="<<ith<<", "<<ID2IRow(ic)<<", "<<ih<<", "<<fHitInputIDs[row.FirstHit()+ih]<<std::endl;     
      if( *fNOutTrackHits>fNHitsTotal ){
	std::cout<<"fNOutTrackHits>fNHitsTotal"<<std::endl;
	//exit(0);
	return;
      }
      out.NHits()++;      
    }    
    if( out.NHits() >= 2 ){
      (*fNOutTracks)++;
    }else {
      (*fNOutTrackHits) = nOutTrackHitsOld;
    }
  }
  timer.Stop();
  fTimers[5]+=timer.CpuTime();
}

GPUh() void AliHLTTPCCATracker::FitTrackFull( AliHLTTPCCATrack &/**/, Float_t * /**/ ) const 
{  
  // fit track with material
#ifdef XXX    
  //* Fit the track   
  FitTrack( iTrack, tt0 );
  if( iTrack.NHits()<=3 ) return;
    
  AliHLTTPCCATrackParam &t = iTrack.Param();
  AliHLTTPCCATrackParam t0 = t;

  t.Chi2() = 0;
  t.NDF() = -5;	
  Bool_t first = 1;

  Int_t iID = iTrack.FirstHitID();
  for( Int_t ih=0; ih<iTrack.NHits(); ih++, iID++ ){
    Int_t *ic = &(fTrackHits[iID]);
    Int_t iRow = ID2IRow(*ic);
    AliHLTTPCCARow &row = fRows[iRow];      
    if( !t0.TransportToX( row.X() ) ) continue;      
    Float_t dy, dz;
    AliHLTTPCCAHit &h = ID2Hit(*ic);

    // check for wrong hits	
    if(0){
      dy = t0.GetY() - h.Y();
      dz = t0.GetZ() - h.Z();
      
      //if( dy*dy > 3.5*3.5*(/*t0.GetErr2Y() + */h.ErrY()*h.ErrY() ) ) continue;//SG!!!
      //if( dz*dz > 3.5*3.5*(/*t0.GetErr2Z() + */h.ErrZ()*h.ErrZ() ) ) continue;      
    }

    if( !t.TransportToX( row.X() ) ) continue;  

    //* Update the track
	    
    if( first ){
      t.Cov()[ 0] = .5*.5;
      t.Cov()[ 1] = 0;
      t.Cov()[ 2] = .5*.5;
      t.Cov()[ 3] = 0;
      t.Cov()[ 4] = 0;
      t.Cov()[ 5] = .2*.2;
      t.Cov()[ 6] = 0;
      t.Cov()[ 7] = 0;
      t.Cov()[ 8] = 0;
      t.Cov()[ 9] = .2*.2;
      t.Cov()[10] = 0;
      t.Cov()[11] = 0;
      t.Cov()[12] = 0;
      t.Cov()[13] = 0;
      t.Cov()[14] = .2*.2;
      t.Chi2() = 0;
      t.NDF() = -5;
    }
    Float_t err2Y, err2Z;
    GetErrors2( iRow, t, err2Y, err2Z );

    if( !t.Filter2( h.Y(), h.Z(), err2Y, err2Z ) ) continue;

    first = 0;      
  }  
  /*
  Float_t cosPhi = iTrack.Param().GetCosPhi();
  p0.Param().TransportToX(ID2Row( iTrack.PointID()[0] ).X());
  p2.Param().TransportToX(ID2Row( iTrack.PointID()[1] ).X());  
  if( p0.Param().GetCosPhi()*cosPhi<0 ){ // change direction
  Float_t *par = p0.Param().Par();
  Float_t *cov = p0.Param().Cov();
  par[2] = -par[2]; // sin phi
  par[3] = -par[3]; // DzDs
    par[4] = -par[4]; // kappa
    cov[3] = -cov[3];
    cov[4] = -cov[4];
    cov[6] = -cov[6];
    cov[7] = -cov[7];
    cov[10] = -cov[10];
    cov[11] = -cov[11];
    p0.Param().CosPhi() = -p0.Param().GetCosPhi();
  }
  */
#endif
}

GPUh() void AliHLTTPCCATracker::FitTrack( AliHLTTPCCATrack &/*track*/, Float_t */*t0[]*/ ) const 
{      
  //* Fit the track   
#ifdef XXX
  AliHLTTPCCAEndPoint &p2 = ID2Point(track.PointID()[1]);
  AliHLTTPCCAHit &c0 = ID2Hit(fTrackHits[p0.TrackHitID()].HitID());	
  AliHLTTPCCAHit &c1 = ID2Hit(fTrackHits[track.HitID()[1]].HitID());	
  AliHLTTPCCAHit &c2 = ID2Hit(fTrackHits[p2.TrackHitID()].HitID());	
  AliHLTTPCCARow &row0 = ID2Row(fTrackHits[p0.TrackHitID()].HitID());
  AliHLTTPCCARow &row1 = ID2Row(fTrackHits[track.HitID()[1]].HitID());
  AliHLTTPCCARow &row2 = ID2Row(fTrackHits[p2.TrackHitID()].HitID());
  Float_t sp0[5] = {row0.X(), c0.Y(), c0.Z(), c0.ErrY(), c0.ErrZ() };
  Float_t sp1[5] = {row1.X(), c1.Y(), c1.Z(), c1.ErrY(), c1.ErrZ() };
  Float_t sp2[5] = {row2.X(), c2.Y(), c2.Z(), c2.ErrY(), c2.ErrZ() };
  //cout<<"Fit track, points ="<<sp0[0]<<" "<<sp0[1]<<" / "<<sp1[0]<<" "<<sp1[1]<<" / "<<sp2[0]<<" "<<sp2[1]<<std::endl;
  if( track.NHits()>=3 ){    
    p0.Param().ConstructXYZ3(sp0,sp1,sp2,p0.Param().CosPhi(), t0);
    p2.Param().ConstructXYZ3(sp2,sp1,sp0,p2.Param().CosPhi(), t0);
    //p2.Param() = p0.Param();
    //p2.Param().TransportToX(row2.X());
    //p2.Param().Par()[1] = -p2.Param().Par()[1];
    //p2.Param().Par()[4] = -p2.Param().Par()[4];
  } else {
    p0.Param().X() = row0.X();
    p0.Param().Y() = c0.Y();
    p0.Param().Z() = c0.Z();
    p0.Param().Err2Y() = c0.ErrY()*c0.ErrY();
    p0.Param().Err2Z() = c0.ErrZ()*c0.ErrZ();
    p2.Param().X() = row2.X();
    p2.Param().Y() = c2.Y();
    p2.Param().Z() = c2.Z();
    p2.Param().Err2Y() = c2.ErrY()*c2.ErrY();
    p2.Param().Err2Z() = c2.ErrZ()*c2.ErrZ();
  }
#endif
}



GPUd() void AliHLTTPCCATracker::GetErrors2( Int_t iRow, const AliHLTTPCCATrackParam &t, Float_t &Err2Y, Float_t &Err2Z ) const
{
  //
  // Use calibrated cluster error from OCDB
  //

  Float_t z = CAMath::Abs((250.-0.275)-CAMath::Abs(t.GetZ()));
  Int_t    type = (iRow<63) ? 0: (iRow>126) ? 1:2;
  Float_t cosPhiInv = CAMath::Abs(t.GetCosPhi())>1.e-2 ?1./t.GetCosPhi() :0;
  Float_t angleY = t.GetSinPhi()*cosPhiInv ;
  Float_t angleZ = t.GetDzDs()*cosPhiInv ;

  Err2Y = fParam.GetClusterError2(0,type, z,angleY);  
  Err2Z = fParam.GetClusterError2(1,type, z,angleZ);
}


#if !defined(HLTCA_GPUCODE)  

GPUh() void AliHLTTPCCATracker::WriteEvent( std::ostream &out ) 
{
  // write event to the file
  for( Int_t iRow=0; iRow<fParam.NRows(); iRow++ ){
    out<<fRows[iRow].FirstHit()<<" "<<fRows[iRow].NHits()<<std::endl;
  } 
  out<<fNHitsTotal<<std::endl;

  Float_t y[fNHitsTotal], z[fNHitsTotal];
 
  for( Int_t iRow=0; iRow<fParam.NRows(); iRow++){
    AliHLTTPCCARow &row = Rows()[iRow];
    Float_t y0 = row.Grid().YMin();
    Float_t z0 = row.Grid().ZMin();
    Float_t stepY = row.HstepY();
    Float_t stepZ = row.HstepZ();
    uint4* tmpint4 = RowData() + row.FullOffset();
    ushort2 *hits = reinterpret_cast<ushort2*>(tmpint4);
    for( Int_t ih=0; ih<fRows[iRow].NHits(); ih++ ){
      Int_t ihTot = row.FirstHit() + ih;
      Int_t id = fHitInputIDs[ihTot];
      ushort2 hh = hits[ih];
      y[id] = y0 + hh.x*stepY;
      z[id] = z0 + hh.y*stepZ;
    }
  }
  for( Int_t ih=0; ih<fNHitsTotal; ih++ ){
    out<<y[ih]<<" "<<z[ih]<<std::endl;
  }
}

GPUh() void AliHLTTPCCATracker::ReadEvent( std::istream &in ) 
{
  //* Read event from file 

  Int_t rowFirstHit[Param().NRows()], rowNHits[Param().NRows()];  
  for( Int_t iRow=0; iRow<Param().NRows(); iRow++ ){
    in>>rowFirstHit[iRow]>>rowNHits[iRow];
  }
  Int_t nHits;
  in >> nHits;
  Float_t y[nHits], z[nHits];
  for( Int_t ih=0; ih<nHits; ih++ ){
    in>>y[ih]>>z[ih];
  }
  ReadEvent( rowFirstHit, rowNHits, y, z, nHits );
} 

GPUh() void AliHLTTPCCATracker::WriteTracks( std::ostream &out ) 
{
  //* Write tracks to file 
 
  out<<fTimers[0]<<std::endl;
  out<<*fNOutTrackHits<<std::endl;
  for( Int_t ih=0; ih<*fNOutTrackHits; ih++ ){
    out<< fOutTrackHits[ih]<<" ";
  }
  out<<std::endl;
  
  out<<*fNOutTracks<<std::endl;

  for( Int_t itr=0; itr<*fNOutTracks; itr++ ){
    AliHLTTPCCAOutTrack &t = fOutTracks[itr];    
    AliHLTTPCCATrackParam &p1 = t.StartPoint();	
    AliHLTTPCCATrackParam &p2 = t.EndPoint();	
    out<< t.NHits()<<" ";
    out<< t.FirstHitRef()<<" ";
    out<< t.OrigTrackID()<<" ";
    out<<std::endl;
    out<< p1.X()<<" ";
    out<< p1.CosPhi()<<" ";
    out<< p1.Chi2()<<" ";
    out<< p1.NDF()<<std::endl;
    for( Int_t i=0; i<5; i++ ) out<<p1.Par()[i]<<" ";
    out<<std::endl;
    for( Int_t i=0; i<15; i++ ) out<<p1.Cov()[i]<<" ";
    out<<std::endl;
    out<< p2.X()<<" ";
    out<< p2.CosPhi()<<" ";
    out<< p2.Chi2()<<" ";
    out<< p2.NDF()<<std::endl;
    for( Int_t i=0; i<5; i++ ) out<<p2.Par()[i]<<" ";
    out<<std::endl;
    for( Int_t i=0; i<15; i++ ) out<<p2.Cov()[i]<<" ";
    out<<std::endl;
  }
}

GPUh() void AliHLTTPCCATracker::ReadTracks( std::istream &in )
{
  //* Read tracks  from file 
  in>>fTimers[0];
  in>>*fNOutTrackHits;  

  for( Int_t ih=0; ih<*fNOutTrackHits; ih++ ){
    in>>fOutTrackHits[ih];
  }
  in>>*fNOutTracks;

  for( Int_t itr=0; itr<*fNOutTracks; itr++ ){
    AliHLTTPCCAOutTrack &t = fOutTracks[itr];    
    AliHLTTPCCATrackParam &p1 = t.StartPoint();	
    AliHLTTPCCATrackParam &p2 = t.EndPoint();	
    in>> t.NHits();
    in>> t.FirstHitRef();
    in>> t.OrigTrackID();    
    in>> p1.X();
    in>> p1.CosPhi();
    in>> p1.Chi2();
    in>> p1.NDF();
    for( Int_t i=0; i<5; i++ ) in>>p1.Par()[i];
    for( Int_t i=0; i<15; i++ ) in>>p1.Cov()[i];
    in>> p2.X();
    in>> p2.CosPhi();
    in>> p2.Chi2();
    in>> p2.NDF();
    for( Int_t i=0; i<5; i++ ) in>>p2.Par()[i];
    for( Int_t i=0; i<15; i++ ) in>>p2.Cov()[i];
  }
}
#endif
