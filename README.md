# FD_Plate

**by Matthew Hamilton (2017)**

C++ Finite Difference Kirschoff Thin Plate Model

***

## Contents
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [FD_Plate](#fdplate)
	- [Contents](#contents)
	- [About](#about)
	- [Build](#build)
		- [Editing](#editing)
		- [Audio Playback](#audio-playback)
	- [Notes](#notes)

<!-- /TOC -->

***

## About

This code illustrates a very basic finite difference plate scheme. No frills, just a very quick implementation for demonstration purposes. Check out the [`FDPlateClass`](https://github.com/mhamilt/FDPlateClass) repository for a more robust implementation.

## Build

A very basic make file is provided. If you open it and read it you will quickly find that make files are a little new to me.

So, you may want to import the `fd-plate.cpp` into your IDE of choice.

### Editing

Edit values in-between the 'Edit Here' banners to change the sound of the plate.

### Audio Playback

At the moment the program simply writes a `.wav` file to a specified directory. The default directory is set to `~/Downloads` (Sorry windows) if now directory is given.

[`cpp-cli-audio-tools`](https://github.com/mhamilt/cpp-cli-audio-tools/tree/bab960fd1d893ff9ba7a92582f1fbd866153b8e0) has been added as a submodule for cli playback in the future.

## Notes

All code tested with Xcode Version 8.2.1 (8C1002) 4.2.1 Compatible Apple LLVM 8.0.0 (clang-800.0.42.1)
