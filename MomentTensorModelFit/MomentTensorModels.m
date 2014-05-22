//
//  MomentTensorModels.m
//  MomentTensorModelFit
//
//  Created by Jeffrey J. Early on 5/22/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import "MomentTensorModels.h"

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