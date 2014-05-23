//
//  DrifterTracksFileReader.h
//  MomentTensorModelFit
//
//  Created by Jeffrey J. Early on 5/23/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface DrifterTracksFileReader : NSObject

- (id) initWithURL: (NSURL *) anURL equation: (GLEquation *) anEquation;

@property(strong, readwrite) NSURL *file;

// Time, in seconds
@property(strong, readwrite) GLFunction *t;

// x position, in meters
@property(strong, readwrite) GLFunction *x;

// y position, in meters
@property(strong, readwrite) GLFunction *y;

@end
