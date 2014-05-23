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
		
		NSArray *result = ellipseComponentsFromMatrixComponents(self.Mxx, self.Myy, self.Mxy);
		//NSArray *result = tensorCompsToEllipseComps(@[self.Mxx, self.Myy, self.Mxy]);
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
		
		NSArray *result = ellipseComponentsFromMatrixComponents(self.Mxx, self.Myy, self.Mxy);
		//NSArray *result = tensorCompsToEllipseComps(@[self.Mxx, self.Myy, self.Mxy]);
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
	}
	return self;
}

- (NSArray *) bestFitToDiffusivityModel
{
	GLFloat kappaScale = 1.0;
	GLScalar *kappa = [GLScalar scalarWithValue: log(1./kappaScale) forEquation: self.equation];
	GLScalar *kappaDelta = [GLScalar scalarWithValue: 2.0 forEquation: self.equation];
		
	GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa] withDeltas: @[kappaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		NSArray *tensorComps = diffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled);
		//NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
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
	GLScalar *sigmaDelta = [GLScalar scalarWithValue: 0.5 forEquation: self.equation];
	
	GLScalar *theta = [GLScalar scalarWithValue: 0.0*M_PI/180. forEquation: self.equation];
	GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: self.equation];
	
	GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, sigma, theta] withDeltas: @[kappaDelta, sigmaDelta, thetaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sigmaUnscaled = [[xArray[1] exponentiate] times: @(sigmaScale)];
		NSArray *tensorComps = strainDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2]);
		NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
		
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
		
	GLScalar *minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
	GLScalar *minSigma = [[minimizer.result[1] exponentiate] times: @(sigmaScale)];
	GLScalar *minTheta = minimizer.result[2];
	GLScalar *minError = minimizer.result[3];
	
	return @[minError, minKappa, minSigma, minTheta];
}

- (NSArray *) bestFitToVorticityStrainDiffusivityModel
{
	GLFloat kappaScale = 0.1;
	GLScalar *kappa = [GLScalar scalarWithValue: log(.1/kappaScale) forEquation: self.equation];
	GLScalar *kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: self.equation];
	
	GLFloat sScale = 1E-6;
	GLScalar *s = [GLScalar scalarWithValue: log(1e-6/sScale) forEquation: self.equation];
	GLScalar *sDelta = [GLScalar scalarWithValue: 0.5 forEquation: self.equation];
	
	GLScalar *theta = [GLScalar scalarWithValue: 0.0*M_PI/180. forEquation: self.equation];
	GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: self.equation];
	
	GLScalar *alpha = [GLScalar scalarWithValue: 0 forEquation: self.equation];
	GLScalar *alphaDelta = [GLScalar scalarWithValue: 1.0 forEquation: self.equation];
	
	// initializing with the results from the previous calculation. we can use log searches if we look for both positive and negative vorticity.
	GLMinimizationOperation *minimizer_positive = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
		GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
		NSArray *tensorComps = strainVorticityDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
		
		EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[self.a, self.b, self.theta]];
		
		return error.result[0];
	}];
	NSArray *posResults = minimizer_positive.result;
	
	alphaDelta = [GLScalar scalarWithValue: -1.0 forEquation: self.equation];
	GLMinimizationOperation *minimizer_negative = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
		GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
		GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
		GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
		GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
		NSArray *tensorComps = strainVorticityDiffusivityModel( self.Mxx0, self.Myy0, self.Mxy0, self.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
		NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
		
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
};

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
	GLScalar * C = [[[zeta_s negate] times: [Mxx1 plus: Myy1]] plus: [[sig_s times: @2] times: Mxy1]];
	
	GLFunction *Maa = [[[exp_s_t times: one_plus_sigma_over_s_div2] times: A] plus: [[exp_minus_s_t times: one_minus_sigma_over_s_div2] times: B]];
	Maa = [Maa plus: [[C times: zeta_s] times: @(0.5)]];
	Maa = [Maa minus: [[[[zeta times: zeta] times: t] plus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
	GLFunction *Mbb = [[[[exp_s_t times: one_minus_sigma_over_s_div2] times: A] plus: [[exp_minus_s_t times: one_plus_sigma_over_s_div2] times: B]] negate];
	Mbb = [Mbb plus: [[C times: zeta_s] times: @(0.5)]];
	Mbb = [Mbb minus: [[[[zeta times: zeta] times: t] minus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
	GLFunction *Mab = [[[[exp_s_t times: zeta_s] times: A] minus: [[exp_minus_s_t times: zeta_s] times: B]] times: @(0.5)];
	Mab = [Mab plus: [[C times: @(0.5)] times: sig_s]];
	Mab = [Mab minus: [[[zeta times: sigma] times: t] times: [[kappa times: @2] dividedBy: s2]]];
	
	GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin times: @(-2.)]]];
	GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin times: @(2.)]]];
	GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
	
	return @[Mxx, Myy, Mxy];
};

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
	
	return @[[lambda_big sqrt], [lambda_small sqrt], theta];
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