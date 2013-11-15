//
//  MomentEllipseFileReader.h
//  MomentTensorModelFit
//
//  Created by Jeffrey Early on 11/15/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface MomentEllipseFileReader : NSObject

- (id) initWithURL: (NSURL *) anURL equation: (GLEquation *) anEquation;

@property(strong, readwrite) NSURL *file;

// Time, in seconds
@property(strong, readwrite) GLFunction *t;

// Semi-major axis, in meters
@property(strong, readwrite) GLFunction *a;

// Semi-minor axis, in meters
@property(strong, readwrite) GLFunction *b;

// Angle of the semi-major axis, in radians
@property(strong, readwrite) GLFunction *angle;

// Initial conditions in terms of the moment tensor components.
@property GLFloat Mxx0;
@property GLFloat Myy0;
@property GLFloat Mxy0;
@end
