//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTTPCGMBORDERTRACK_H
#define ALIHLTTPCGMBORDERTRACK_H

#include "AliGPUCommonDefGPU.h"

/**
 * @class AliGPUTPCGMBorderTrack
 *
 * The class describes TPC slice tracks at sector borders.
 * Used in AliGPUTPCGMMerger
 *
 */
class AliGPUTPCGMBorderTrack
{

  public:
	struct Range
	{
		int fId;
		float fMin, fMax;
		static bool CompMin(const Range &a, const Range &b) { return a.fMin < b.fMin; }
		static bool CompMax(const Range &a, const Range &b) { return a.fMax < b.fMax; }
	};

	GPUd() int TrackID() const { return fTrackID; }
	GPUd() short NClusters() const { return fNClusters; }
	GPUd() short Row() const { return fRow; }
	GPUd() const float *Par() const { return fP; }
	GPUd() float ZOffset() const { return fZOffset; }
	GPUd() const float *Cov() const { return fC; }
	GPUd() const float *CovD() const { return fD; }

	GPUd() void SetTrackID(int v) { fTrackID = v; }
	GPUd() void SetNClusters(short v) { fNClusters = v; }
	GPUd() void SetRow(short v) { fRow = v; }
	GPUd() void SetPar(int i, float x) { fP[i] = x; }
	GPUd() void SetZOffset(float v) { fZOffset = v; }
	GPUd() void SetCov(int i, float x) { fC[i] = x; }
	GPUd() void SetCovD(int i, float x) { fD[i] = x; }

	GPUd() static bool CheckChi2(float x1, float y1, float cx1, float cxy1, float cy1,
	                      float x2, float y2, float cx2, float cxy2, float cy2, float chi2cut)
	{
		//* Calculate Chi2/ndf deviation
		float dx = x1 - x2;
		float dy = y1 - y2;
		float cx = cx1 + cx2;
		float cxy = cxy1 + cxy2;
		float cy = cy1 + cy2;
		float det = cx * cy - cxy * cxy;
		//printf("Res %f Det %f Cut %f %s   -   ", ( cy*dx - (cxy+cxy)*dy )*dx + cx*dy*dy, det, (det + det) * chi2cut, (CAMath::Abs(( cy*dx - (cxy+cxy)*dy )*dx + cx*dy*dy) < CAMath::Abs((det+det)*chi2cut)) ? "OK" : "Fail");
		return (CAMath::Abs((cy * dx - (cxy + cxy) * dy) * dx + cx * dy * dy) < CAMath::Abs((det + det) * chi2cut));
	}

	GPUd() bool CheckChi2Y(const AliGPUTPCGMBorderTrack &t, float chi2cut) const
	{
		float d = fP[0] - t.fP[0];
		return (d * d < chi2cut * (fC[0] + t.fC[0]));
	}

	GPUd() bool CheckChi2Z(const AliGPUTPCGMBorderTrack &t, float chi2cut) const
	{
		float d = fP[1] - t.fP[1] + (fZOffset - t.fZOffset);
		return (d * d < chi2cut * (fC[1] + t.fC[1]));
	}

	GPUd() bool CheckChi2QPt(const AliGPUTPCGMBorderTrack &t, float chi2cut) const
	{
		float d = fP[4] - t.fP[4];
		return (d * d < chi2cut * (fC[4] + t.fC[4]));
	}

	GPUd() bool CheckChi2YS(const AliGPUTPCGMBorderTrack &t, float chi2cut) const
	{
		return CheckChi2(fP[0], fP[2], fC[0], fD[0], fC[2],
		                 t.fP[0], t.fP[2], t.fC[0], t.fD[0], t.fC[2], chi2cut);
	}

	GPUd() bool CheckChi2ZT(const AliGPUTPCGMBorderTrack &t, float chi2cut) const
	{
		return CheckChi2(fP[1], fP[3], fC[1], fD[1], fC[3],
		                 t.fP[1] + (t.fZOffset - fZOffset), t.fP[3], t.fC[1], t.fD[1], t.fC[3], chi2cut);
	}

  private:
	int fTrackID;     // track index
	short fNClusters; // n clusters
	short fRow;
	float fP[5];
	float fZOffset;
	float fC[5];
	float fD[2];
};

#endif
