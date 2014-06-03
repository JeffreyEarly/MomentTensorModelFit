//
//  EllipseErrorOperation.m
//  MomentTensorModelFit
//
//  Created by Jeffrey Early on 11/15/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import "EllipseErrorOperation.h"
#include "ellipse_ellipse_overlap.h"

#define ERROR_METHOD 1

@implementation EllipseErrorOperation

- (EllipseErrorOperation *) initWithParametersFromEllipseA: (NSArray *) paramsA ellipseB: (NSArray *) paramsB
{
    GLEquation *equation = [paramsA[0] equation];
    GLDimension *tDim = [paramsA[0] dimensions][0];
    NSUInteger nPoints = tDim.nPoints;
    
    NSArray *result = @[[GLScalar scalarWithType: kGLRealDataFormat forEquation: equation]];
    
    NSMutableArray *operand = [NSMutableArray arrayWithArray: paramsA];
    [operand addObjectsFromArray: paramsB];
    
    NSMutableArray *buffers = [NSMutableArray array];
    buffers[0] = [[GLBuffer alloc] initWithLength: nPoints*sizeof(GLFloat)];
    buffers[1] = [[GLBuffer alloc] initWithLength: nPoints*sizeof(GLFloat)];
    
    variableOperation errorOperation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        GLFloat *a_obs = (GLFloat *) [operandArray[0] bytes];
        GLFloat *b_obs = (GLFloat *) [operandArray[1] bytes];
        GLFloat *angle_obs = (GLFloat *) [operandArray[2] bytes];
        GLFloat *a_model = (GLFloat *) [operandArray[3] bytes];
        GLFloat *b_model = (GLFloat *) [operandArray[4] bytes];
        GLFloat *angle_model = (GLFloat *) [operandArray[5] bytes];
        
        GLFloat *totalError = (GLFloat *) [resultArray[0] bytes];
        
        GLFloat *errorAtTime = (GLFloat *) [bufferArray[0] bytes];
        GLFloat *areaAtTime = (GLFloat *) [bufferArray[1] bytes];
		
#if ERROR_METHOD == 0 || ERROR_METHOD == 1
        
        dispatch_apply(tDim.nPoints, dispatch_get_global_queue (DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^(size_t iTime) {
			double X[4], Y[4];
			int nroots, rtnCode;
			double overlapArea = ellipse_ellipse_overlap(angle_obs[iTime], a_obs[iTime], b_obs[iTime], 0.0, 0.0, angle_model[iTime], a_model[iTime], b_model[iTime],0.0, 0.0,X, Y, &nroots, &rtnCode);
			double obsArea =  M_PI * a_obs[iTime] *  b_obs[iTime];
			double modelArea = M_PI * a_model[iTime] *  b_model[iTime];
			errorAtTime[iTime] = (modelArea - overlapArea) + (obsArea - overlapArea);
			areaAtTime[iTime] = obsArea;
        });
		// There are several logical ways to weight the error.
#if ERROR_METHOD == 0
        // 1. sum (error[i]/area[i])
		// For Brownian motion this is sub-optimal and results in over-weighting the first ellipses.
        vGL_vdiv( areaAtTime, 1, errorAtTime, 1, errorAtTime, 1, nPoints);
		errorAtTime[0] = 0.0;
        vGL_sve( errorAtTime, 1, totalError, nPoints);
        *totalError = (*totalError)/nPoints;
#else
		// 2. (sum error[i])/(sum area[i])
        // This has the advantage that the total error is weighted towards the larger ellipse at the end---good for Brownian motion
        GLFloat errorSum;
        vGL_sve( errorAtTime, 1, &errorSum, nPoints);

        GLFloat areaSum;
        vGL_sve( areaAtTime, 1, &areaSum, nPoints);

        *totalError = errorSum/areaSum;
#endif

#elif ERROR_METHOD == 2
		size_t iTime = tDim.nPoints -1;
		double X[4], Y[4];
		int nroots, rtnCode;
		double overlapArea = ellipse_ellipse_overlap(angle_obs[iTime], a_obs[iTime], b_obs[iTime], 0.0, 0.0, angle_model[iTime], a_model[iTime], b_model[iTime],0.0, 0.0,X, Y, &nroots, &rtnCode);
		double obsArea =  M_PI * a_obs[iTime] *  b_obs[iTime];
		double modelArea = M_PI * a_model[iTime] *  b_model[iTime];
		errorAtTime[iTime] = (modelArea - overlapArea) + (obsArea - overlapArea);
		areaAtTime[iTime] = obsArea;
		
		*totalError = errorAtTime[iTime]/areaAtTime[iTime];
		
#endif
    };
    
    if ((self=[super initWithResult: result operand: operand buffers: buffers operation:errorOperation])) {
        
    }
    
    return self;
}

@end
