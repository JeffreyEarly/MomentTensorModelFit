//
//  MomentTensorModels.h
//  MomentTensorModelFit
//
//  Created by Jeffrey J. Early on 5/22/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface MomentTensorModels : NSObject

- (MomentTensorModels *) initWithXPositions: (GLFunction *) x yPositions: (GLFunction *) y time: (GLFunction *) t;

- (MomentTensorModels *) initWithMxx: (GLFunction *) Mxx Myy: (GLFunction *) Myy Mxy: (GLFunction *) Mxy time: (GLFunction *) t;

- (MomentTensorModels *) initWithA: (GLFunction *) a b: (GLFunction *) b theta: (GLFunction *) theta time: (GLFunction *) t;

// Returns error, kappa
- (NSArray *) bestFitToDiffusivityModel;

// Returns error, kappa, zeta
- (NSArray *) bestFitToVorticityDiffusivityModel;

// Returns error, kappa, sigma, theta
- (NSArray *) bestFitToStrainDiffusivityModel;

// Returns error, kappa, sigma, theta, zeta
- (NSArray *) bestFitToVorticityStrainDiffusivityModel;
- (NSArray *) bestFitToVorticityStrainDiffusivityModelWithStartPoint: (NSArray *) startPoint;

- (NSArray *) bestFitToVorticityStrainMatchedDiffusivityModelWithStartPoint: (NSArray *) startPoint;
- (NSArray *) bestFitToVorticityStrainDominatedDiffusivityModelWithStartPoint: (NSArray *) startPoint;
- (NSArray *) bestFitToVorticityDominatedStrainDiffusivityModelWithStartPoint: (NSArray *) startPoint;
@end



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

NSArray * vorticityDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *zeta);
// This model assume that you're handing it sigma and zeta of the same magnitude (zeta can be negative, however).
NSArray * strainVorticityMatchedDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta, GLScalar *zeta );

NSArray * strainVorticityDominantedDiffusivityModel(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta, GLScalar *zeta );

// Simple little utility function that converts tensor components, into ellipse components
NSArray *tensorCompsToEllipseComps( NSArray *tensorComp );

// Returns a, b, theta---semi-major, semi-minor, and angle of semi-major axis.
NSArray *ellipseComponentsFromMatrixComponents( GLFunction *Mxx, GLFunction *Myy, GLFunction *Mxy);