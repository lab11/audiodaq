#!/usr/bin/python
import sys
import os
import time

import audiotools as at
import struct
from collections import deque

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class audioDaqProcessor:
    def __init__(self):
        # The deque container has O(1) inserting and removal.

        # Buffer holds 300 rolling points, deque
        # was chosen because it is a lot more efficient
        # to pop items in and out of it than a list.
        self.buff = deque()

        # Difference between each sequential set of points
        # zero is added as the first value to prevent errors
        # when the algorithm is bootstrapping itself
        self.diff = deque()
        self.diff.append(0)

        # We keep track of the current state in
        # self.state, and use it to reference a function
        # pointer to the function being called. Super simple
        # way to process data.
        self.stateMachine = {
            "START": self.fnStart,
            "E1": self.fnE1,
            "E2": self.fnE2,
            "E3": self.fnE3
        }

        self.state = "START"


        self.frameMaxCounts = [0, 0, 0, 0]
        self.edgeCount = 0

        self.pointCount = 0

        self.maxDiff = 0
        self.minDiff = 0

        self.diffDetectThres1 = 0.75
        self.diffDetectThres2 = 0.65
        self.diffDetectThres1 = 0.70
        self.diffDetectThres2 = 0.65
        self.window = 10

        self.upPeak1 = 300
        self.downPeak1 = 300
        self.upPeak2 = 300
        self.downPeak2 = 300

        # contains the sequential sample index of the
        # detected and inferred edges.
        self.points = open('points.csv', 'w')

        # contains the four estimated values per frame
        # stored with their locations
        self.values = open('values.csv', 'w')

        self.envolope = open('envelope.csv', 'w')

        # Output data stored in this file in the format
        # time (in seconds), output value (in volts)
        self.out = open('out.csv', 'w')

    # public function called by the object that has
    # instanced this class, pushes a single audio sample
    # into the algorithm for processing.
    def receive(self, pcmValue):
        self.pointCount += 1
        self.buff.append(pcmValue)
        self.updateDiff()
        if len(self.buff) > 300:
            self.buff.popleft()
            self.process()

        # Debug print
        self.envolope.write('%d,%d,%d,%d\n' % (pcmValue, self.diff[-1], self.maxDiff, self.minDiff))

    # public function to close the debug and output file handles
    def close(self):
        self.points.close()
        self.envolope.close()
        self.out.close()
        self.values.close()

    # the location of the different audio peaks is
    # stored as a location in the rolling buffer, we want
    # them to move with the buffer, so we increment them
    # if this fails, we reset state, obviously we didn't
    # detect a full frame in time
    def movePointers(self):
        self.upPeak1 -= 1
        self.downPeak1 -= 1
        self.upPeak2 -= 1
        self.downPeak2 -= 1
        if self.upPeak1 < 0:
            self.resetState()
            print 'Resetting state'

    def resetState(self):
        self.upPeak1 = 300
        self.downPeak1 = 300
        self.upPeak2 = 300
        self.downPeak2 = 300
        self.state = "START"

    # called internally, calls the outFn function pointer
    # if present with the recorded data
    def sendData(self, t, x):
        if hasattr(self, 'outFn'):
            self.outFn(t, x)

    # core invocation function
    def process(self):
        self.movePointers()
        self.stateMachine[self.state]()

    def fnE1(self):
        # if the rate of rise of the current middle of the buffer
        # is less than a lower bound, we've found the first down
        # peak.
        if self.diff[150] < self.minDiff * self.diffDetectThres2:
            self.downPeak1 = self.getIndexOf(min)
            self.edgeCount += 1
            #self.points.write('%d\n' % (self.pointCount - 150))

            # we now estimate the width of the space between
            # peaks, and the location of the next two peaks, which
            # would be tough to detect due to the unknowable value
            # of the analog input signal
            peakDelta = self.downPeak1 - self.upPeak1

            self.upPeak2 = self.downPeak1 + peakDelta
            self.downPeak2 = self.upPeak2 + peakDelta
            self.state = "E2"

            if peakDelta < 8:
                self.state = "START"
            #Debug print
            print 'delta: %d u-edge: %d d-edge: %d' % (peakDelta, self.upPeak1, self.downPeak1)

    def fnE2(self):
        # once we've shifted in enough data to find the
        # last down peak, we move to a processing stage.
        if self.downPeak2 == 150:
            self.state = "E3"

    def fnE3(self):
        # now we're finding what is considered the last edge of the
        # previous frame, and the first edge of the next frame.

        # if this is true, we've estimated correctly, else we silently
        # fail and let the frame shifting expire and state transition
        # back to start.
        if self.diff[150] > self.maxDiff * self.diffDetectThres1:
            newPeak = self.getIndexOf(max)
            val1 = self.getValueOf(self.upPeak1, self.downPeak1, 'v')
            val2 = self.getValueOf(self.downPeak1, self.upPeak2, 'g')
            val3 = self.getValueOf(self.upPeak2, self.downPeak2, 's')
            val4 = self.getValueOf(self.downPeak2, newPeak, 'g')
            try:
                t = float(self.pointCount) / 44100
                x = float(val3 - float((val2 + val4) / 2)) / float(val1 - ((val2 + val4 / 2))) * 1.8
            except:
                return
            self.sendData(t, x)

            # write the final two points to the detected set
            self.points.write('%d\n' % (self.pointCount - (300 - self.upPeak1)))
            self.points.write('%d\n' % (self.pointCount - (300 - self.downPeak1)))
            self.points.write('%d\n' % (self.pointCount - (300 - self.upPeak2)))
            self.points.write('%d\n' % (self.pointCount - (300 - self.downPeak2)))

            self.resetState()
            self.state = "E1"
            self.upPeak1 = newPeak

            # Debug print
            self.out.write('%f,%f\n' % (t, x))
            print '(%d, %d, %d, %d)' % (val1, val2, val3, val4)

    def fnStart(self):
        # Goal: Find first large peak up, frame data
        print 'maxDiff: %d minDiff: %d' % (self.maxDiff, self.minDiff)
        self.resetState()

        # if the rate of rise of the current middle of the buffer
        # is greater than a threshold, we've found the first up
        # peak.
        if self.diff[150] > self.maxDiff * self.diffDetectThres1:
            self.upPeak1 = self.getIndexOf(max)
            self.state = "E1"
            self.edgeCount += 1

            #Debug print
            print 'e1: %d' % self.pointCount
            #self.points.write('%d\n' % (self.pointCount - 150))

    def getValueOf(self, start, stop, pt):
        if stop - start < 20:
            return 0

        movement = int(round((stop - start) * 0.1))
        skip = int(round((stop - start) * 0.02))

        estVal = 0
        vals = self.getSubList(self.buff, start + skip, stop - skip)

        if pt == 'v':
            estVal = sum(vals) / len(vals)
            estVal = min(vals)
        elif pt == 'g':
            estVal = min(vals)
            estVal = estVal + sum(vals) / len(vals)
            estVal = estVal / 2
        elif pt == 's':
            estVal = sum(vals) / len(vals)

        self.values.write('%f,%d\n' % (float(self.pointCount - 300 + start + skip + movement / 2) / 44100, estVal))
        return estVal

    def getIndexOf(self, fnMinMax):
            sub = self.getSubList(self.diff, 150 - self.window / 2, 150 + self.window / 2)
            iMax = 150 - self.window / 2 + sub.index(fnMinMax(sub))
            return iMax

    def getSubList(self, dQ, start, stop):
        ret = []
        for i in range(start, stop):
            ret.append(dQ[i])
        return ret

    def updateDiff(self):
        buffLen = len(self.buff)
        if buffLen < 2:
            return
        newDiff = self.buff[buffLen - 1] - self.buff[buffLen - 2]

        self.diff.append(self.buff[buffLen - 1] - self.buff[buffLen - 2])

        if self.maxDiff < newDiff: self.maxDiff = newDiff
        if self.minDiff > newDiff: self.minDiff = newDiff

        # This is computationally expensive, so we only do it
        # if there's a good chance of the max/min changing.
        if len(self.diff) > 300:
            old = self.diff[0]
            self.diff.popleft()
            if self.maxDiff == old: self.maxDiff = max(self.diff)
            if self.minDiff == old: self.minDiff = min(self.diff)

# Main code section


def doDecode(fileName):
    print 'Opening input data auio stream... decoding'
    aF = at.open(fileName)
    pcmAf = aF.to_pcm()

    a = audioDaqProcessor()

    xs = []
    ts = []

    def handleData(t, x):
        xs.append(x)
        ts.append(t)

    a.outFn = handleData
    rawFile = open('raw.csv', 'w')
    hackCounter = 0

    while True:
        frame = pcmAf.read(256)
        for i in range(0, frame.frames):
            byteArray = frame.channel(0).frame(i).to_bytes(False, True)
            pcmVal = struct.unpack('h', byteArray)
            rawFile.write(str(float(hackCounter) / float(44100)) + "," + str(pcmVal[0]) + '\n')
            hackCounter = hackCounter + 1
            a.receive(pcmVal[0])
        if frame.frames < 64:
            print 'End of file found. Breaking.'
            break

    a.close()

    rawFile.close()

    out = open('out1.csv', 'w')
    for i in range(0, len(xs)):
        out.write('%f,%f\n' % (ts[i], xs[i]))
    out.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ts, xs)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('voltage (V)')
    ax.set_title('Recovered Signal')
    ax.grid(True)

    fig.savefig('output.png')


while True:
    if os.path.exists('/uploads/file1.m4a'):
        os.system("mv %s %s" % ('/uploads/file1.m4a', 'file1.m4a')
        doDecode('file1.m4a')
        time.sleep(1)

