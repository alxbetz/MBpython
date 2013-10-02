import numpy as np
import math

def get2dAngle(omitAxis,ax):
	if omitAxis == "z":
		tmp = [ax[0],ax[1],0]
	elif omitAxis == "y":
		tmp = [ax[0],0,ax[2]]
	elif omitAxis == "x":
		tmp [0,ax[1],ax[2]]
	return math.acos(np.dot(tmp,[1,0,0]))
		

def mRotX(angle):
	return np.array([[1,0,0,1],[0,math.cos(angle),math.sin(angle),0],[0,-(math.sin(angle)),math.cos(angle),0],[0,0,0,1]])

def mRotY(angle):
	return np.array([[math.cos(angle),0,-(math.sin(angle)),0],[0,1,0,1][math.sin(angle),0,math.cos(angle),0],[0,0,0,1]])

def mRotZ(angle):
	return np.array([math.cos(angle),math.sin(angle),0,0],[-(math.sin(angle)),math.cos(angle),0,0],[0,0,1,0][0,0,0,1]])

def rotateM(axis,angle):
	angleZ = get2dAngle("z",axis)
	angleY = get2dAngle("y",axis)
	res = np.dot(mRotZ(angleZ),mRotY(angleY))
	res = np.dot(res,mRotX(angle))
	res = np.dot(res,mRotY(-angleY))
	res = np.dot(res,mRotZ(-angleZ))
	return res

def rotate(xyz,axis,angle):
	rotM = rotateM(axis,angle)
	for v in xyz:
		v= np.dot(v,rotM)
			

