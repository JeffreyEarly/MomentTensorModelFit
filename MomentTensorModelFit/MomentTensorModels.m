//
//  MomentTensorModels.m
//  MomentTensorModelFit
//
//  Created by Jeffrey J. Early on 5/22/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import "MomentTensorModels.h"
#import <GLNumericalModelingKit/GLOperationOptimizer.h>
#import "EllipseErrorOperation.h"

@interface MomentTensorModels ()
@property GLFunction *Mxx;
@property GLFunction *Myy;
@property GLFunction *Mxy;
@property GLFloat Mxx0;
@property GLFloat Myy0;
@property GLFloat Mxy0;
@property GLFloat MxxN;
@property GLFloat MyyN;
@property GLFloat MxyN;
@property GLFloat maxT;
@property GLFunction *a;
@property GLFunction *b;
@property GLFunction *theta;
@property GLFunction *t;
@property GLEquation *equation;
@end

@implementation MomentTensorModels

- (MomentTensorModels *) initWithXPositions: (GLFunction *) x yPositions: (GLFunction *) y time: (GLFunction *) t
{
	if ((self=[super init])) {
		
		GLFunction *x_com = [x mean: 1];
		GLFunction *y_com = [y mean: 1];
		
		GLFunction *q = [x plus: [x_com negate]];
		GLFunction *r = [y plus: [y_com negate]];
		
		self.Mxx = [[q times: q] mean: 1];
		self.Myy = [[r times: r] mean: 1];
		self.Mxy = [[q times: r] mean: 1];
		self.t = t;
		self.equation = x.equation;
		self.Mxx0 = self.Mxx.pointerValue[0];
		self.Myy0 = self.Myy.pointerValue[0];
		self.Mxy0 = self.Mxy.pointerValue[0];
        
        NSUInteger n=self.t.nDataPoints-1;
        self.MxxN = self.Mxx.pointerValue[n];
        self.MyyN = self.Myy.pointerValue[n];
        self.MxyN = self.Mxy.pointerValue[n];
        self.maxT = self.t.pointerValue[n];
		
//		self.Mxx = [[self.Mxx minus: @(self.Mxx0)] abs];
//		self.Myy = [[self.Myy minus: @(self.Myy0)] abs];
//		self.Mxy = [[self.Mxy minus: @(self.Mxy0)] abs];
//		
//		GLFunction *q0 = [GLFunction functionOfRealTypeWithDimensions: @[q.dimensions[1]] forEquation:self.equation];
//		GLFunction *r0 = [GLFunction functionOfRealTypeWithDimensions: @[r.dimensions[1]] forEquation:self.equation];
//		for (NSUInteger i=0; i<q0.nDataPoints; i++) {
//			q0.pointerValue[i] = q.pointerValue[i];
//			r0.pointerValue[i] = r.pointerValue[i];
//		}
//
//		GLFunction *q_bar = [q plus: [q0 negate]];
//		GLFunction *r_bar = [r plus: [r0 negate]];
//		self.Mxx = [[q_bar times: q_bar] mean: 1];
//		self.Myy = [[r_bar times: r_bar] mean: 1];
//		self.Mxy = [[q_bar times: r_bar] mean: 1];
		
		NSArray *result = ellipseComponentsFromMatrixComponents(self.Mxx, self.Myy, self.Mxy);
		self.a = result[0];
		self.b = result[1];
		self.theta = result[2];
	}
	return self;
}

- (MomentTensorModels *) initWithMxx: (GLFunction *) Mxx Myy: (GLFunction *) Myy Mxy: (GLFunction *) Mxy time: (GLFunction *) t
{
	if ((self=[super init])) {
		self.Mxx = Mxx;
		self.Myy = Myy;
		self.Mxy = Mxy;
		self.Mxx0 = self.Mxx.pointerValue[0];
		self.Myy0 = self.Myy.pointerValue[0];
		self.Mxy0 = self.Mxy.pointerValue[0];
		self.t = t;
		self.equation = self.Mxx.equation;
		
        NSUInteger n=self.t.nDataPoints-1;
        self.MxxN = self.Mxx.pointerValue[n];
        self.MyyN = self.Myy.pointerValue[n];
        self.MxyN = self.Mxy.pointerValue[n];
        self.maxT = self.t.pointerValue[n];
        
		NSArray *result = ellipseComponentsFromMatrixComponents(self.Mxx, self.Myy, self.Mxy);
		self.a = result[0];
		self.b = result[1];
		self.theta = result[2];
	}
	return self;
}

- (MomentTensorModels *) initWithA: (GLFunction *) a b: (GLFunction *) b theta: (GLFunction *) theta time: (GLFunction *) t
{
	if ((self=[super init])) {
		self.a = a;
		self.b = b;
		self.theta = theta;
		self.t = t;
		self.equation = self.a.equation;
		
		GLFloat D2p = pow(self.a.pointerValue[0],2.0);
        GLFloat D2m = pow(self.b.pointerValue[0],2.0);
        GLFloat alpha = self.theta.pointerValue[0];
        
        self.Mxx0 = D2p*cos(alpha)*cos(alpha) + D2m*sin(alpha)*sin(alpha);
        self.Myy0 = D2p*sin(alpha)*sin(alpha) + D2m*cos(alpha)*cos(alpha);
        self.Mxy0 = (D2p-D2m)*sin(alpha)*cos(alpha);
        
        NSUInteger n=self.t.nDataPoints-1;
        D2p = pow(self.a.pointerValue[n],2.0);
        D2m = pow(self.b.pointerValue[n],2.0);
        alpha = self.theta.pointerValue[n];
        self.MxxN = D2p*cos(alpha)*cos(alpha) + D2m*sin(alpha)*sin(alpha);
        self.MyyN = D2p*sin(alpha)*sin(alpha) + D2m*cos(alpha)*cos(alpha);
        self.MxyN = (D2p-D2m)*sin(alpha)*cos(alpha);
        self.maxT = self.t.pointerValue[n];
	}
	return self;
}

- (NSArray *) bestFitToDiffusivityModel
{
	GLFloat kappaScale = 0.1;
	GLScalar *kappa = [GLScalar scalarWithValue: log(0.1/kappaScale) forEquation: self.equation];
	GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
	
	GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa] withDeltas: @[kappaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		NSArray *tensorComps = diffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
		
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	
	GLScalar *minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
	GLScalar *minError = minimizer.result[1];
	
	return @[minError, minKappa];
}

- (NSArray *) bestFitToStrainDiffusivityModel
{
	// kappaDelta is carefully chosen. It needs to represent a sort of 'size of parameter space' that we want to explore.
	// So here, we make it move around in fairly big chunks, like orders of magnitude.
	// Yup, big steps are the best, for ALL parameters.
	GLFloat kappaScale = 0.1;
	GLScalar *kappa = [GLScalar scalarWithValue: log(.1/kappaScale) forEquation: self.equation];
	GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
	
	GLFloat sigmaScale = 4.0E-6;
	GLScalar *sigma = [GLScalar scalarWithValue: log(1E-6/sigmaScale) forEquation: self.equation];
	GLScalar *sigmaDelta = [GLScalar scalarWithValue: 1.0 forEquation: self.equation];
	
	GLScalar *theta = [GLScalar scalarWithValue: 0.0*M_PI/180. forEquation: self.equation];
	GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: self.equation];
	
	GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, sigma, theta] withDeltas: @[kappaDelta, sigmaDelta, thetaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sigmaUnscaled = [[xArray[1] exponentiate] times: @(sigmaScale)];
		NSArray *tensorComps = strainDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2]);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
		
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
		
	GLScalar *minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
	GLScalar *minSigma = [[minimizer.result[1] exponentiate] times: @(sigmaScale)];
	GLScalar *minTheta = minimizer.result[2];
	GLScalar *minError = minimizer.result[3];
	
	return @[minError, minKappa, minSigma, minTheta];
}

- (NSArray *) bestFitToStrainDiffusivityModelWithFixedStrainAngle: (GLScalar *) theta0
{
    // kappaDelta is carefully chosen. It needs to represent a sort of 'size of parameter space' that we want to explore.
    // So here, we make it move around in fairly big chunks, like orders of magnitude.
    // Yup, big steps are the best, for ALL parameters.
    GLFloat kappaScale = 0.1;
    GLScalar *kappa = [GLScalar scalarWithValue: log(.1/kappaScale) forEquation: self.equation];
    GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
    
    GLFloat sigmaScale = 4.0E-6;
    GLScalar *sigma = [GLScalar scalarWithValue: log(1E-6/sigmaScale) forEquation: self.equation];
    GLScalar *sigmaDelta = [GLScalar scalarWithValue: 1.0 forEquation: self.equation];
    
    GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, sigma] withDeltas: @[kappaDelta, sigmaDelta] forFunction:^(NSArray *xArray) {
        GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
        GLScalar *sigmaUnscaled = [[xArray[1] exponentiate] times: @(sigmaScale)];
        NSArray *tensorComps = strainDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, theta0);
        NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
        
        EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
        
        return error.result[0];
    }];
    
    GLScalar *minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
    GLScalar *minSigma = [[minimizer.result[1] exponentiate] times: @(sigmaScale)];
    GLScalar *minTheta = theta0;
    GLScalar *minError = minimizer.result[2];
    
    return @[minError, minKappa, minSigma, minTheta];
}

- (NSArray *) bestFitToVorticityDiffusivityModel
{
    // We use the first and last moments to estimate a good set of starting parameters
    GLFloat kappaStart = ((self.MxxN+self.MyyN)-(self.Mxx0+self.Myy0))/(4*self.maxT);
    kappaStart = kappaStart < 0 ? 1e-4 : kappaStart;
    
    GLFloat B = -2.0*self.Mxy0;
    GLFloat C = self.Mxx0 - self.Myy0;
    
    GLFloat arg = ((self.MxxN-self.MyyN)/C) + (2*self.MxyN/B);
    arg *= B*C/(B*B+C*C);
    if (arg < -1) {
        arg = -1;
    } else if (arg > 1) {
        arg = 1;
    }
    GLFloat zetaStart = fabs(asin(arg)/self.maxT);
    
    GLFloat kappaScale = 0.1;
    GLScalar *kappa = [GLScalar scalarWithValue: log(kappaStart/kappaScale) forEquation: self.equation];
    GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
    
    GLFloat zetaScale = 0.1/self.maxT;
    GLScalar *zeta = [GLScalar scalarWithValue: log(zetaStart/zetaScale) forEquation: self.equation];
    GLScalar *zetaDelta = [GLScalar scalarWithValue: 1.0 forEquation: self.equation];
    
    // First check the positive values of zeta
    GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, zeta] withDeltas: @[kappaDelta, zetaDelta] forFunction:^(NSArray *xArray) {
        GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
        GLScalar *zetaUnscaled = [[xArray[1] exponentiate] times: @(zetaScale)];
        NSArray *tensorComps = vorticityDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, zetaUnscaled);
        NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
        
        EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
        
        return error.result[0];
    }];
    
    GLScalar *minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
    GLScalar *minZeta = [[minimizer.result[1] exponentiate] times: @(zetaScale)];
    GLScalar *minError = minimizer.result[2];
    
    // Now check the negative values of zeta
    minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, zeta] withDeltas: @[kappaDelta, zetaDelta] forFunction:^(NSArray *xArray) {
        GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
        GLScalar *zetaUnscaled = [[xArray[1] exponentiate] times: @(-zetaScale)];
        NSArray *tensorComps = vorticityDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, zetaUnscaled);
        NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
        
        EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
        
        return error.result[0];
    }];
    
    GLScalar *minKappaNeg = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
    GLScalar *minZetaNeg = [[minimizer.result[1] exponentiate] times: @(-zetaScale)];
    GLScalar *minErrorNeg = minimizer.result[2];
    
    if (minError.pointerValue[0] < minErrorNeg.pointerValue[0]) {
        return @[minError, minKappa, minZeta];
    } else {
        return @[minErrorNeg, minKappaNeg, minZetaNeg];
    }
}

- (NSArray *) bestFitToVorticityStrainDiffusivityModel
{
    GLScalar *kappa = [GLScalar scalarWithValue: 0.1 forEquation: self.equation];
    GLScalar *sScale = [GLScalar scalarWithValue: 1e-5 forEquation: self.equation];
    GLScalar *theta = [GLScalar scalarWithValue: 0.0 forEquation: self.equation];
    return [self bestFitToVorticityStrainDiffusivityModelWithStartPoint: @[kappa, sScale, theta]];
}

- (NSArray *) bestFitToVorticityStrainDiffusivityModelWithStartPoint: (NSArray *) startPoint
{
	NSArray * strainDominant = [self bestFitToVorticityStrainDominatedDiffusivityModelWithStartPoint: startPoint];
	NSArray * matched = [self bestFitToVorticityStrainMatchedDiffusivityModelWithStartPoint: startPoint];
	NSArray * vorticityDominant = [self bestFitToVorticityDominatedStrainDiffusivityModelWithStartPoint: startPoint];
	
	NSArray *results = @[strainDominant, matched,vorticityDominant];
	NSArray *sortedResults = [results sortedArrayUsingComparator: ^NSComparisonResult(NSArray *obj1, NSArray *obj2) {
		GLScalar *error1 = obj1[0];
		GLScalar *error2 = obj2[0];
		if (error1.pointerValue[0] < error2.pointerValue[0]) {
			return NSOrderedAscending;
		} else {
			return NSOrderedDescending;
		}
	}];
	
	return sortedResults[0];
}

- (NSArray *) bestFitToVorticityStrainDominatedDiffusivityModelWithStartPoint: (NSArray *) startPoint
{
	GLFloat kappaScale = 0.1;
	GLScalar *kappa = [[startPoint[0] times: @(1./kappaScale)] log];;
	GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
	
	GLFloat sScale = 1E-5;
	GLScalar *s = [[startPoint[1] times: @(1./sScale)] log];
	GLScalar *sDelta = [GLScalar scalarWithValue: 1.0 forEquation: self.equation];
	
	GLScalar *theta = startPoint[2];
	GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: self.equation];
	
	GLScalar *alpha = [GLScalar scalarWithValue: 0.0 forEquation: self.equation];
	GLScalar *alphaDelta = [GLScalar scalarWithValue: 1 forEquation: self.equation];
	
	// initializing with the results from the previous calculation. we can use log searches if we look for both positive and negative vorticity.
	GLMinimizationOperation *minimizer_positive = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
		GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
		NSArray *tensorComps = strainVorticityDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
		
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	NSArray *posResults = minimizer_positive.result;
	
	alpha = [GLScalar scalarWithValue: 0.0 forEquation: self.equation];
	alphaDelta = [GLScalar scalarWithValue: -1 forEquation: self.equation];
	GLMinimizationOperation *minimizer_negative = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
		GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
		NSArray *tensorComps = strainVorticityDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );

		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	NSArray *negResults = minimizer_negative.result;
	
	GLScalar *minPos = posResults[4];
	GLScalar *minNeg = negResults[4];
	GLFloat pos = *(minPos.pointerValue);
	GLFloat neg = *(minNeg.pointerValue);
	
	NSArray *results = pos < neg ? posResults : negResults;
	
	GLScalar *minKappa = [[results[0] exponentiate] times: @(kappaScale)];
	GLScalar *minS = [[results[1] exponentiate] times: @(sScale)];
	GLScalar *minSigma = [minS times: [results[3] cosh]];
	GLScalar *minTheta = results[2];
	GLScalar *minZeta = [minS times: [results[3] sinh]];
	GLScalar *minError = results[4];
	
	return @[minError, minKappa, minSigma, minTheta, minZeta];
}

- (NSArray *) bestFitToVorticityStrainMatchedDiffusivityModelWithStartPoint: (NSArray *) startPoint
{
	GLFloat kappaScale = 0.1;
	GLScalar *kappa = [[startPoint[0] times: @(1./kappaScale)] log];;
	GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
	
	GLFloat sScale = 1E-5;
	GLScalar *s = [[startPoint[1] times: @(1./sScale)] log];
	GLScalar *sDelta = [GLScalar scalarWithValue: 1.0 forEquation: self.equation];
	
	GLScalar *theta = startPoint[2];
	GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: self.equation];
	
	GLMinimizationOperation *minimizer_special_pos = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta] withDeltas: @[kappaDelta, sDelta, thetaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sigmaUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *zetaUnscaled = sigmaUnscaled;
		NSArray *tensorComps = strainVorticityMatchedDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	NSArray *specPosResults = minimizer_special_pos.result;
	
	GLMinimizationOperation *minimizer_special_neg = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta] withDeltas: @[kappaDelta, sDelta, thetaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sigmaUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *zetaUnscaled = [sigmaUnscaled negate];
		NSArray *tensorComps = strainVorticityMatchedDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	NSArray *specNegResults = minimizer_special_neg.result;
	
	GLScalar *minSpecPos = specPosResults[3];
	GLScalar *minSpecNeg = specNegResults[3];
	GLFloat specPos = *(minSpecPos.pointerValue);
	GLFloat specNeg = *(minSpecNeg.pointerValue);
	
	NSArray *results = specPos < specNeg ? specPosResults : specNegResults;
	
	GLScalar *minKappa = [[results[0] exponentiate] times: @(kappaScale)];
	GLScalar *minSigma = [[results[1] exponentiate] times: @(sScale)];
	GLScalar *minTheta = results[2];
	GLScalar *minZeta = specPos < specNeg ? minSigma : [minSigma negate];
	GLScalar *minError = results[3];
	
	return @[minError, minKappa, minSigma, minTheta, minZeta];
}

- (NSArray *) bestFitToVorticityDominatedStrainDiffusivityModelWithStartPoint: (NSArray *) startPoint
{
	GLFloat kappaScale = 0.1;
	GLScalar *kappa = [[startPoint[0] times: @(1./kappaScale)] log];;
	GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
	
	// Really \bar{s}
	GLFloat sScale = 1E-5;
	GLScalar *s = [[startPoint[1] times: @(1./sScale)] log];
	GLScalar *sDelta = [GLScalar scalarWithValue: 1.0 forEquation: self.equation];
	
	GLScalar *theta = startPoint[2];
	GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: self.equation];
	
	GLScalar *alpha = [GLScalar scalarWithValue: 0 forEquation: self.equation];
	GLScalar *alphaDelta = [GLScalar scalarWithValue: 0.5 forEquation: self.equation];
	
	// initializing with the results from the previous calculation. we can use log searches if we look for both positive and negative vorticity.
	GLMinimizationOperation *minimizer_positive = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *sigmaUnscaled = [[sUnscaled times: [xArray[3] sinh]] abs];
		GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] cosh]];
		NSArray *tensorComps = strainVorticityDominantedDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
		
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	NSArray *posResults = minimizer_positive.result;
	
//	alpha = [GLScalar scalarWithValue: 0 forEquation: self.equation];
//	alphaDelta = [GLScalar scalarWithValue: 0.5 forEquation: self.equation];
	GLMinimizationOperation *minimizer_negative = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *sigmaUnscaled = [[sUnscaled times: [xArray[3] sinh]] abs];
		GLScalar *zetaUnscaled = [[sUnscaled times: [xArray[3] cosh]] negate];
		NSArray *tensorComps = strainVorticityDominantedDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = ellipseComponentsFromMatrixComponents( tensorComps[0], tensorComps[1], tensorComps[2] );
		
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	NSArray *negResults = minimizer_negative.result;
	
	GLScalar *minPos = posResults[4];
	GLScalar *minNeg = negResults[4];
	GLFloat pos = *(minPos.pointerValue);
	GLFloat neg = *(minNeg.pointerValue);
	
	NSArray *results = pos < neg ? posResults : negResults;
	
	GLScalar *minKappa = [[results[0] exponentiate] times: @(kappaScale)];
	GLScalar *minS = [[results[1] exponentiate] times: @(sScale)];
	GLScalar *minSigma = [[minS times: [results[3] sinh]] abs];
	GLScalar *minTheta = results[2];
	GLScalar *minZeta = pos < neg ? [minS times: [results[3] cosh]] : [[minS times: [results[3] cosh]] negate];
	GLScalar *minError = results[4];
	
	return @[minError, minKappa, minSigma, minTheta, minZeta];
}

@end

// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
// It takes one parameter, kappa.
// It returns (Mxx, Myy, Mxy) for all time t.
NSArray *diffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa )
{
	GLFunction *Mxx = [[t times: [kappa times: @(2)]] plus: @(Mxx0)];
	GLFunction *Myy = [[t times: [kappa times: @(2)]] plus: @(Myy0)];
	GLFunction *Mxy = [[t times: @(0)] plus: @(Mxy0)];
	
	return @[Mxx, Myy, Mxy];
}

// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
// It takes three parameters, kappa, sigma, and theta.
// It returns (Mxx, Myy, Mxy) for all time t.
NSArray * vorticityDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *zeta)
{
	GLFloat A = 0.5*(Mxx0 + Myy0);
	GLFloat B = -Mxy0;
	GLFloat C = 0.5*(Mxx0 - Myy0);
	
	GLFunction *sin_zeta_t = [[t times: zeta] sin];
	GLFunction *cos_zeta_t = [[t times: zeta] cos];
	
	GLFunction * tks = [t times:[kappa times: @(2)]];
	
	GLFunction *Mxx = [[[tks plus: [sin_zeta_t times: @(B)]] plus: [cos_zeta_t times: @(C)]] plus: @(A)];
	GLFunction *Myy = [[[tks minus: [sin_zeta_t times: @(B)]] minus: [cos_zeta_t times: @(C)]] plus: @(A)];
	GLFunction *Mxy = [[sin_zeta_t times: @(C)] plus: [cos_zeta_t times:@(-B)]];
	
	return @[Mxx, Myy, Mxy];
}

// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
// It takes three parameters, kappa, sigma, and theta.
// It returns (Mxx, Myy, Mxy) for all time t.
NSArray * strainDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta )
{
	GLScalar *cos_t = [theta cos];
	GLScalar *sin_t = [theta sin];
	
	GLScalar *cos2 = [cos_t times: cos_t];
	GLScalar *sin2 = [sin_t times: sin_t];
	GLScalar *cossin = [cos_t times: sin_t];
	
	GLScalar * tks = [[kappa times: @(2)] dividedBy: sigma];
	
	GLScalar * A = [[[cos2 times: @(Mxx0)] plus: [sin2 times: @(Myy0)]] plus: [[cossin times: @(2.*Mxy0)] plus: tks]];
	GLScalar * B = [[[sin2 times: @(Mxx0)] plus: [cos2 times: @(Myy0)]] minus: [[cossin times: @(2.*Mxy0)] plus: tks]];
	GLScalar * C = [[[cossin times: @(-Mxx0)] plus: [cossin times: @(Myy0)]] plus: [[cos2 minus: sin2] times: @(Mxy0)]];
	
	GLFunction *Maa = [[[[t times: sigma] exponentiate] times: A] minus: tks];
	GLFunction *Mbb = [[[[t times: [sigma negate]] exponentiate] times: B] plus: tks];
	GLFunction *Mab = [[t scalarMultiply: 0.0] plus: C];
	
	GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin times: @(-2.)]]];
	GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin times: @(2.)]]];
	GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
	
	return @[Mxx, Myy, Mxy];
};

// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
// It takes four parameters, kappa, sigma, theta, and zeta.
// It returns (Mxx, Myy, Mxy) for all time t.
NSArray * strainVorticityDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta, GLScalar *zeta )
{
	GLScalar *cos_t = [theta cos];
	GLScalar *sin_t = [theta sin];
	
	GLScalar *cos2 = [cos_t times: cos_t];
	GLScalar *sin2 = [sin_t times: sin_t];
	GLScalar *cossin = [cos_t times: sin_t];
	
	GLScalar * s2 = [[sigma times: sigma] minus: [zeta times: zeta]];
	GLScalar * s = [s2 sqrt];
	GLScalar * tks = [[[kappa times: @(2)] times: sigma] dividedBy: s2];
	GLScalar * sig_s = [sigma dividedBy: s];
	GLScalar * zeta_s = [zeta dividedBy: s];
	
	GLScalar * one_plus_sigma_over_s_div2 = [[sig_s plus: @(1)] times: @(0.5)];
	GLScalar * one_minus_sigma_over_s_div2 = [[sig_s minus: @(1)] times: @(-0.5)];
	GLFunction * exp_s_t = [[t times: s] exponentiate];
	GLFunction * exp_minus_s_t = [[t times: [s negate]] exponentiate];
	
	GLScalar * Mxx1 = [[[cos2 times: @(Mxx0)] plus: [sin2 times: @(Myy0)]] plus: [cossin times: @(2.*Mxy0)]];
	GLScalar * Myy1 = [[[sin2 times: @(Mxx0)] plus: [cos2 times: @(Myy0)]] minus: [cossin times: @(2.*Mxy0)]];
	GLScalar * Mxy1 = [[[cossin times: @(-Mxx0)] plus: [cossin times: @(Myy0)]] plus: [[cos2 minus: sin2] times: @(Mxy0)]];
	
	GLScalar * A = [[[[one_plus_sigma_over_s_div2 times: Mxx1] minus: [zeta_s times: Mxy1]] minus: [one_minus_sigma_over_s_div2 times: Myy1]] plus: tks];
	GLScalar * B = [[[[one_minus_sigma_over_s_div2 times: Mxx1] plus: [zeta_s times: Mxy1]] minus: [one_plus_sigma_over_s_div2 times: Myy1]] plus: tks];
	GLScalar * C = [[[[zeta negate] times: [Mxx1 plus: Myy1]] plus: [[sigma times: @2] times: Mxy1]] dividedBy: s2];
	
	GLFunction *Maa = [[[exp_s_t times: one_plus_sigma_over_s_div2] times: A] plus: [[exp_minus_s_t times: one_minus_sigma_over_s_div2] times: B]];
	Maa = [Maa plus: [[C times: zeta] times: @(0.5)]];
	Maa = [Maa minus: [[[[zeta times: zeta] times: t] plus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
	GLFunction *Mbb = [[[[exp_s_t times: one_minus_sigma_over_s_div2] times: A] plus: [[exp_minus_s_t times: one_plus_sigma_over_s_div2] times: B]] negate];
	Mbb = [Mbb plus: [[C times: zeta] times: @(0.5)]];
	Mbb = [Mbb minus: [[[[zeta times: zeta] times: t] minus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
	GLFunction *Mab = [[[[exp_s_t times: zeta_s] times: A] minus: [[exp_minus_s_t times: zeta_s] times: B]] times: @(0.5)];
	Mab = [Mab plus: [[C times: @(0.5)] times: sigma]];
	Mab = [Mab minus: [[[zeta times: sigma] times: t] times: [[kappa times: @2] dividedBy: s2]]];
	
	GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin times: @(-2.)]]];
	GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin times: @(2.)]]];
	GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
	
	return @[Mxx, Myy, Mxy];
};

// Special case of the above where zeta^2 = sigma^2
NSArray * strainVorticityMatchedDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta, GLScalar *zeta )
{
	GLScalar *cos_t = [theta cos];
	GLScalar *sin_t = [theta sin];
	
	GLScalar *cos2 = [cos_t times: cos_t];
	GLScalar *sin2 = [sin_t times: sin_t];
	GLScalar *cossin = [cos_t times: sin_t];
	
	GLScalar * Mxx1 = [[[cos2 times: @(Mxx0)] plus: [sin2 times: @(Myy0)]] plus: [cossin times: @(2.*Mxy0)]];
	GLScalar * Myy1 = [[[sin2 times: @(Mxx0)] plus: [cos2 times: @(Myy0)]] minus: [cossin times: @(2.*Mxy0)]];
	GLScalar * Mxy1 = [[[cossin times: @(-Mxx0)] plus: [cossin times: @(Myy0)]] plus: [[cos2 minus: sin2] times: @(Mxy0)]];
	
	GLScalar * A = [[sigma times:[Mxx1 plus: Myy1]] minus: [[zeta times: @(2)] times: Mxy1]];
	GLScalar * B = [Mxx1 minus: Myy1];
	GLScalar * C = [Mxx1 plus: Myy1];
	
	GLFunction *ksigmatcubed = [[[[t times: t] times:t] times: [[sigma times: sigma] times: kappa] ] times: @(1./3.)];
	GLFunction *ksigmatsquared = [[t times: t] times: [sigma times: kappa]];
	GLFunction *twokappat = [[t times: kappa] times: @(2)];
	GLFunction *sigmatsquared = [[t times: t] times: sigma];
	
	GLFunction *Maa = [[ksigmatcubed plus: ksigmatsquared] plus: twokappat];
	Maa = [Maa plus: [[A times: [sigmatsquared plus: [t times: @(2)]]] times: @(0.25)]];
	Maa = [Maa plus: [[B times: [[t times: sigma] plus: @(1)]] times: @(0.5)]];
	Maa = [Maa plus: [C times: @(0.5)]];
	
	GLFunction *Mbb = [[ksigmatcubed minus: ksigmatsquared] plus: twokappat];
	Mbb = [Mbb plus: [[A times: [sigmatsquared minus: [t times: @(2)]]] times: @(0.25)]];
	Mbb = [Mbb plus: [[B times: [[t times: sigma] minus: @(1)]] times: @(0.5)]];
	Mbb = [Mbb plus: [C times: @(0.5)]];
	
	GLFunction *Mab = ksigmatcubed;
	GLScalar *half = [GLScalar scalarWithValue: 0.5 forEquation: Mab.equation];
	Mab = [Mab plus: [A times: [[[zeta times: [t times: t]] times: @(0.25)] minus: [half dividedBy: zeta]]]];
	Mab = [Mab plus: [B times: [[t times: zeta] times: @(0.5)]]];
	Mab = [Mab plus: [C times: [[sigma dividedBy: zeta] times: @(0.5)]]];
	
	GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin times: @(-2.)]]];
	GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin times: @(2.)]]];
	GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
	
	return @[Mxx, Myy, Mxy];
}

// The case where vorticity is stronger than strain
NSArray * strainVorticityDominantedDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta, GLScalar *zeta )
{
	GLScalar *cos_t = [theta cos];
	GLScalar *sin_t = [theta sin];
	
	GLScalar *cos2 = [cos_t times: cos_t];
	GLScalar *sin2 = [sin_t times: sin_t];
	GLScalar *cossin = [cos_t times: sin_t];
	
	GLScalar * s2 = [[zeta times: zeta] minus: [sigma times: sigma]];
	GLScalar * s = [s2 sqrt];
	GLScalar * tks = [[[kappa times: @(2)] times: sigma] dividedBy: s2];
	GLScalar * sig_s = [sigma dividedBy: s];
	GLScalar * zeta_s = [zeta dividedBy: s];
	GLScalar * sig_zeta = [sigma dividedBy: zeta];
	
	GLFunction * sin_s_t = [[t times: s] sin];
	GLFunction * cos_s_t = [[t times: s] cos];
	
	GLScalar * Mxx1 = [[[cos2 times: @(Mxx0)] plus: [sin2 times: @(Myy0)]] plus: [cossin times: @(2.*Mxy0)]];
	GLScalar * Myy1 = [[[sin2 times: @(Mxx0)] plus: [cos2 times: @(Myy0)]] minus: [cossin times: @(2.*Mxy0)]];
	GLScalar * Mxy1 = [[[cossin times: @(-Mxx0)] plus: [cossin times: @(Myy0)]] plus: [[cos2 minus: sin2] times: @(Mxy0)]];
	
	GLScalar * Aover2 = [[[Mxx1 minus: Myy1] times: @(0.5)] minus: tks];
	GLScalar * Bover2 = [[[[[Mxx1 plus: Myy1] times: sigma] times: @(0.5)] minus: [zeta times: Mxy1]] dividedBy: s];
	GLScalar * Cover2 = [[[[[Mxx1 plus: Myy1] times: [zeta times: zeta]] times: @(0.5)] minus: [[zeta times: sigma] times: Mxy1]] dividedBy: s2];
	
	GLFunction *Maa = [[cos_s_t plus: [sig_s times: sin_s_t]] times: Aover2];
	Maa = [Maa plus: [[sin_s_t minus: [sig_s times: cos_s_t]] times: Bover2]];
	Maa = [Maa plus: Cover2];
	Maa = [Maa plus: [[[[zeta times: zeta] times: t] plus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
	
	GLFunction *Mbb = [[[cos_s_t minus: [sig_s times: sin_s_t]] times: Aover2] negate];
	Mbb = [Mbb minus: [[sin_s_t plus: [sig_s times: cos_s_t]] times: Bover2]];
	Mbb = [Mbb plus: Cover2];
	Mbb = [Mbb minus: [[[[zeta times: zeta] times: t] minus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
	
	GLFunction *Mab = [[[sin_s_t times: zeta_s] times: Aover2] minus: [[cos_s_t times: zeta_s] times: Bover2]];
	Mab = [Mab plus: [Cover2 times: sig_zeta]];
	Mab = [Mab minus: [[[zeta times: sigma] times: t] times: [[kappa times: @2] dividedBy: s2]]];
	
	GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin times: @(-2.)]]];
	GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin times: @(2.)]]];
	GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
	
	return @[Mxx, Myy, Mxy];
}

NSArray *ellipseComponentsFromMatrixComponents( GLFunction *Mxx, GLFunction *Myy, GLFunction *Mxy)
{
	// For a real symmetric matrix [a b; b d], the eigenvalues are
	// 2*lambda = (a+d) \pm \sqrt{ (a-d)^2 + 4b^2 }
	GLFunction *a_plus_d = [Mxx plus: Myy];
	GLFunction *a_minus_d = [Mxx minus: Myy];
	GLFunction *radical = [[[a_minus_d times: a_minus_d] plus: [[Mxy times: Mxy] times: @4]] sqrt];
	GLFunction *lambda_big = [[a_plus_d plus: radical] times: @(0.5)];
	GLFunction *lambda_small = [[a_plus_d minus: radical] times: @(0.5)];
	
	//GLFunction *theta = [[lambda_big minus: Mxx] atan2: Myy];
	GLFunction *theta = [Mxy atan2: [lambda_big minus: Myy]];
	
	return @[[[lambda_big abs] sqrt], [[lambda_small abs] sqrt], theta];
}

// Simple little utility function that converts tensor components, into ellipse components
NSArray *tensorCompsToEllipseComps( NSArray *tensorComp )
{
	GLFunction *Mxx = tensorComp[0];
	GLFunction *Myy = tensorComp[1];
	GLFunction *Mxy = tensorComp[2];
	
	GLFunction *descriminant = [[[[Mxx minus: Myy] times: [Mxx minus: Myy]] plus: [[Mxy scalarMultiply: 2] times:[Mxy scalarMultiply: 2]]] sqrt];
	GLFunction *D2_plus = [[[Mxx plus: Myy] plus: descriminant] scalarMultiply: 0.5];
	GLFunction *D2_minus = [[[Mxx plus: Myy] minus: descriminant] scalarMultiply: 0.5];
	
	GLFunction *semiMajor = [D2_plus sqrt];
	GLFunction *semiMinor = [D2_minus sqrt];
	GLFunction *angle = [[[D2_plus minus: Mxx] dividedBy: Mxy] atan];
	
	return @[semiMajor, semiMinor, angle];
};
