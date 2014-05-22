//
//  MomentTensorModels.h
//  MomentTensorModelFit
//
//  Created by Jeffrey J. Early on 5/22/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>


// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
// It takes one parameter, kappa.
// It returns (Mxx, Myy, Mxy) for all time t.
NSArray *diffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa );

// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
// It takes three parameters, kappa, sigma, and theta.
// It returns (Mxx, Myy, Mxy) for all time t.
NSArray * strainDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta );

// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
// It takes four parameters, kappa, sigma, theta, and zeta.
// It returns (Mxx, Myy, Mxy) for all time t.
NSArray * strainVorticityDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta, GLScalar *zeta );

// Simple little utility function that converts tensor components, into ellipse components
NSArray *tensorCompsToEllipseComps( NSArray *tensorComp );