import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.core.shape_base import vstack
from sys import exit

#Initialisieren der Landmarks
LM1 = [3.5,  2]
LM2 = [3.5, -2]
LM3 = [  0, -4]
LMR = 0.5 #Radius der Landmarks

'''
 Definiert die Flächen eine Snapshot/ einer Retina-Abbildung (Circle) als 2 Punkte 
 und die Rotationswinkel im Bogenmaß zur Mitte beider Punkte.
 Der Mittelpunkt des Kreises wird bei der Berechnung des 
 Mittelpunktes der Fläche verwendet.
'''
class Area(object):
    def __init__(self, p1, p2, circlecenter):
        self.p1 = p1
        self.p2 = p2

        theta = polarConverter2(p1, p2, circlecenter)
        self.center = (theta[0] + ((theta[1]-theta[0])/2)) % (2*math.pi)
    def __repr__(self) -> str:
        return str(self.p1) + " | " + str(self.p2) + " | " + str(self.center)
    def __eq__(self, o: object) -> bool:
        return self.p1 == object.p1 and self.p2 == object.p2 and self.center == object.center


'''
 Definiert den Kreis eines Snapshots/einer Retina-Abbildung
 als eine Liste von Flächen, gegenübergesetzten Flächen und
 der Mittepunkt des Kreises 
'''
class Circle(object):
    def __init__(self, areas, center):
        self.center = center
        self.areas = []
        self.op_areas = []
        # Sort Areas
        theta = []
        order = list(range(len(areas)))
        for a in areas:
            theta.append( polarConverter(a, center)[0])
        for t in theta:
            for i in theta:
                if t > i:
                    temp = t
                    t = i
                    i = temp
        for o in order:
            self.areas.append(areas[o])
        # Calculate oposite areas
        if len(self.areas) == 3:
            self.op_areas.append(Area(self.areas[0].p2, self.areas[2].p1, center))
            self.op_areas.append(Area(self.areas[2].p2, self.areas[1].p1, center))
            self.op_areas.append(Area(self.areas[1].p2, self.areas[0].p1, center))
        if len(self.areas) == 2:
            self.op_areas.append(Area(self.areas[0].p1, self.areas[1].p2, center))
            self.op_areas.append(Area(self.areas[1].p1, self.areas[0].p2, center))
        if len(self.areas) == 1:
            self.op_areas.append(Area(self.areas[0].p2, self.areas[0].p1, center))

'''
 Mathematische Funktionen zur Bearbeitung von Vektoren (hier: orthogonale Rotation und Längenberechnung)
'''
# Vector Rotation
def rotateVector2DClockwise(vector):
    return np.array([vector[1], -vector[0]])
def rotateVector2DCounterClockwise(vector):
    return np.array([-vector[1], vector[0]])
def vector2DLenght(vector):
    return math.sqrt(vector[0]**2 + vector[1]**2)

#Skalarprodukt
def skalar(u, v):
    return u[0]*v[0]+u[1]*v[1]

'''
 Gibt die Fläche zurück, die durch ein Landmark auf einen Kreis mit gegebenem Radius 
 projiziert wird. Dessen Mittelpunkt liegt auf middle.
'''
def calculateArea(lm, middle, radius):
    #Convert to numpy Arrays (Vectors)
    lm = np.array(lm)
    middle = np.array(middle)
    M_LM = lm - middle # vector mitte zu LM
    M_LMR = rotateVector2DClockwise(M_LM) #vector mitte zu LM rotiert
    M_LMR_05 = M_LMR / (2 * vector2DLenght(M_LMR)) # Vektor durch Länge teilen, dann anpassen (mit *2) wegen r=0.5 von LM
    M_LM_EDGE0 = M_LM + M_LMR_05 # raender von LM1 
    M_LM_EDGE1 = M_LM - M_LMR_05 # andere ^
    M_LM_EDGE0_P1 = (M_LM_EDGE0 / vector2DLenght(M_LM_EDGE0)) * radius #center zu schnittpunkt Kreis Rand1
    M_LM_EDGE1_P2 = (M_LM_EDGE1 / vector2DLenght(M_LM_EDGE1)) * radius # ^ andere
    return Area(middle + M_LM_EDGE0_P1, middle +M_LM_EDGE1_P2, middle)

'''
 Beide Funktionen geben die Rotation im Bogenmass ’Theta' der beiden Punkte
 die einen Fläche bilden zurück. Verwendet wird der Mittelpunkt des Kreises middle für die Berechnung.
'''
def polarConverter(area, middle):
    theta = [0,0]
    theta[0] = math.atan2(area.p1[1] - middle[1], area.p1[0] - middle[0])
    theta[1] = math.atan2(area.p2[1] - middle[1], area.p2[0] - middle[0]) 

    if theta[0] < 0:
        theta[0] += 2*(math.pi)
    if theta[1] < 0:
        theta[1] += 2*(math.pi)    
    if theta[0] > theta[1]:
        theta[1] += 2*(math.pi)

    return theta
def polarConverter2(p1, p2, middle):
    theta = [0,0]
    theta[0] = math.atan2(p1[1] - middle[1], p1[0] - middle[0])
    theta[1] = math.atan2(p2[1] - middle[1], p2[0] - middle[0]) 

    if theta[0] < 0:
        theta[0] += 2*(math.pi) 
    if theta[1] < 0:
        theta[1] += 2*(math.pi)   
    if theta[0] > theta[1]:
        theta[1] += 2*(math.pi)

    return theta
def polarConverter3(punkt, middle):
    return math.atan2(punkt[1]-middle[1], punkt[0]-middle[0])
'''
 Funktion zur Bestimmung, ob sich zwei gegebene Flächen überlappen.
 middle ist der Mittelpunkt des Kreises, zu dem die Flächen gehören.
'''
def areaOverlap(area1, area2, middle):
    # Polar Coordinates
    angles = [[0,0],[0,0]]
    angles[0] = polarConverter(area1, middle)
    angles[1] = polarConverter(area2, middle)

    # a1.p1 in a2
    a = (angles[0][0] >= angles[1][0] and angles[0][0] <= angles[1][1]) 
    # a1.p2 in a2
    b = (angles[0][1] >= angles[1][0] and angles[0][1] <= angles[1][1]) 
    # a2.p1 in a1
    c = (angles[1][0] >= angles[0][0] and angles[1][0] <= angles[0][1])
    # a2.p2 in a1
    d = (angles[1][1] >= angles[0][0] and angles[1][1] <= angles[0][1])
    
    return a or b or c or d

'''
 Gibt einen Circle zurück, der eine Retina-Abbildung oder einen Snapshot repräsentiert, mit einem Mittelpunkt middle und einem gegebenen Radius radius.
 Bei der Berechnung der entsprechenden Flächen wird geprüft, ob sich Flächen überlappen und verbindet sie automatisch, wenn dies der Fall ist. 
'''
def takeSnapshot(middle, radius):
    areas = [] # Real Projected areas
    valid_areas = [] # Joined reduced areas
    
    #Calculate Areas
    areas.append(calculateArea(LM1, middle, radius))
    areas.append(calculateArea(LM2, middle, radius))
    areas.append(calculateArea(LM3, middle, radius))
    
    # Calculate the rotation degrees in radians of the points of all the areas projected by the landmarks
    theta = [0 for x in range(6)] 
    for i in range(0, 3):
            theta[i*2] = math.atan2(areas[i].p1[1]- middle[1], areas[i].p1[0]- middle[0])
            theta[i*2 + 1] = math.atan2(areas[i].p2[1]- middle[1], areas[i].p2[0]- middle[0])
            
            if theta[i*2] < 0:
                theta[i*2] += 2*(math.pi)
            
            if theta[i*2 + 1] < 0:
                theta[i*2 + 1] += 2*(math.pi)
                
            if theta[i*2] > theta[i*2 + 1]:
                theta[i*2 + 1] += 2*(math.pi)

    area1 = True
    area2 = True
    area3 = True
    # Überprüfen mit areaOverlap ob es Überlappungen gibt und verwenden von theta, um zu bestimmen, wie die verbunden werden sollen.
    if (areaOverlap(areas[0], areas[1], middle) and areaOverlap(areas[1], areas[2], middle)) or (areaOverlap(areas[0], areas[2], middle) and areaOverlap(areas[2], areas[1], middle)):
        area1 = False
        area2 = False
        area3 = False
        index_min = 0
        index_max = 0
        for x in theta:
            if x == min(theta):
                break
            index_min += 1
        for x in theta:
            if x == max(theta):
                break
            index_max += 1
        valid_areas.append(Area(areas[int(index_min/2)].p1, areas[int(index_max/2)].p2, middle))
    else:
        if areaOverlap(areas[0], areas[1], middle):
            area1 = False
            area2 = False
            # A1 innerhalb von a A2
            if (theta[0] >= theta[2] and theta[0] <= theta[3]) and (theta[1] >= theta[2] and theta[1] <= theta[3]):
                valid_areas.append(areas[1])
            # A2 innerhalb von a A1
            elif (theta[2] >= theta[0] and theta[2] <= theta[1]) and (theta[3] >= theta[0] and theta[3] <= theta[1]):
                valid_areas.append(areas[0])
            # A1 rechts von A2 a1.p1 zb. a2 | a2.p1 -> a1.p2
            elif theta[0] >= theta[2] and theta[0] <= theta[3]:
                valid_areas.append(Area(areas[1].p1, areas[0].p2, middle))
            # A1 links von A2 a1.p2 zb. a2
            elif theta[1] >= theta[2] and theta[1] <= theta[3]:
                valid_areas.append(Area(areas[0].p1, areas[1].p2, middle))
            #print("Overlap A1-A2")
        if areaOverlap(areas[0], areas[2], middle):
            area1 = False
            area3 = False
            # A1 innerhalb von a A3
            if (theta[0] >= theta[4] and theta[0] <= theta[5]) and (theta[1] >= theta[4] and theta[1] <= theta[5]):
                valid_areas.append(areas[2])
            # A3 innerhalb von a A1
            elif (theta[4] >= theta[0] and theta[4] <= theta[1]) and (theta[5] >= theta[0] and theta[5] <= theta[1]):
                valid_areas.append(areas[0])
            # A1 rechts von A3 a1.p1 zb. a3 | a2.p1 -> a1.p2
            elif theta[0] >= theta[4] and theta[0] <= theta[5]:
                valid_areas.append(Area(areas[2].p1, areas[0].p2, middle))
            # A1 links von A3 a1.p2 zb. a3
            elif theta[1] >= theta[4] and theta[1] <= theta[5]:
                valid_areas.append(Area(areas[0].p1, areas[2].p2, middle))
            #print("Overlap A1-A3")
        if areaOverlap(areas[1], areas[2], middle):
            area2 = False
            area3 = False
            # A3 innerhalb von a A2
            if (theta[4] >= theta[2] and theta[4] <= theta[3]) and (theta[5] >= theta[2] and theta[5] <= theta[3]):
                valid_areas.append(areas[1])
            # A2 innerhalb von a A3
            elif (theta[2] >= theta[4] and theta[2] <= theta[5]) and (theta[3] >= theta[4] and theta[3] <= theta[5]):
                valid_areas.append(areas[2])
            # A3 rechts von A2 a1.p1 zb. a2 | a2.p1 -> a1.p2
            elif theta[4] >= theta[2] and theta[4] <= theta[3]:
                valid_areas.append(Area(areas[1].p1, areas[2].p2, middle))
            # A3 links von A2 a1.p2 zb. a2
            elif theta[5] >= theta[2] and theta[5] <= theta[3]:
                valid_areas.append(Area(areas[2].p1, areas[1].p2, middle))
            #print("Overlap A2-A3")
    if area1:
        valid_areas.append(areas[0])
    if area2:
        valid_areas.append(areas[1])
    if area3:
        valid_areas.append(areas[2])
    
    return Circle(valid_areas, middle) 

# Main
'''
Snapshot wird bei Startpunkt [0, 0] initialisiert
'''
snapshot = takeSnapshot(np.array([0,0]), 1)

'''
Menü mit drei Auswahlmöglichkeiten öffnet sich für den Benutzer.
'''
while True:    
    #Plot wird eine Größe zugeordnet
    fig, ax = plt.subplots(figsize = (7,7))
    
    print("\nBitte wählen Sie einen Option: \n(Erwartete eingaben: '1' '2' '0')")
    print("1. Einen Homing-Vektor ausgeben")
    print("2. Alle Homing-Vektoren ausgeben")
    print("0. Programm beenden\n")
    while True:
        try:
            eingabe = int(input())
        except:
            eingabe = 3
        if eingabe != 1 or eingabe != 2 or eingabe != 0:
            break
        else:
            print("Falsches Eingabe! Bitte versuchen Sie es erneut.")
    if eingabe == 0:
        exit()
    elif eingabe == 1:
        while True:
            print("Nur Koordinaten im bereich [-7;7] sind erlaubt")
            print("Bitte geben Sie Ihre x-Koordinate ein:")
            try: 
                x_coord = int(input())
                print("Bitte geben Sie Ihre y-Koordinate ein:")
                y_coord = int(input())
            except :
                continue

            if ((x_coord >= -7 and x_coord <= 7) and (y_coord >= -7 and y_coord <= 7)) or\
                ((x_coord != 0 and y_coord != 0) or (x_coord != LM1[0] and y_coord != LM1[1]) or\
                (x_coord != LM2[0] and y_coord != LM2[1]) or (x_coord != LM3[0] and y_coord != LM3[1])):
                break
            else:
                print("Falsche Eingabe! Bitte versuchen Sie es erneut.")
        #print('X: ' + str(j-7) + " | Y: " + str(l-7))
                
        '''
        Retina-Abbildung wird bei den gegebenen Koordinaten erstellt.
        '''
        retina = takeSnapshot(np.array([x_coord, y_coord]), 2)
        
        '''
        Rotations- und Translations-Vektoren werden durch Zuordnung und Vergleich der Flächen bestimmt.
        '''
        rotation_vectors = []
        translation_vectors = []
        # Berechnen des Rotations- und Translationsvektor jedes Snapshot-Fläche und seiner gepaarten Retina-Fläche 
        for s_area in snapshot.areas:
            abstand = []
            pair = 0
            #Retina-Abbildungsflächen werden den passenden Snapshot-Flächen zugeordnet
            for r_area in retina.areas:
                berechnung = abs(s_area.center - r_area.center)
                
                if (berechnung > (math.pi)):
                    berechnung = (2*math.pi - berechnung)
                abstand.append(abs(berechnung))
                
                if (abs(berechnung)) == min(abstand):
                    pair = r_area
                    
            # Rotationsvektoren werden erstellt
            if s_area.center < pair.center:
                #print(r_area.center)
                rotation_vectors.append(rotateVector2DClockwise([math.cos(pair.center), math.sin(pair.center)]))
            elif s_area.center > pair.center:
                rotation_vectors.append(rotateVector2DCounterClockwise([math.cos(pair.center), math.sin(pair.center)]))
            else:
                rotation_vectors.append(np.array([0,0]))
                
            # Translationsvektoren werden erstellt
            # Breite der Areas berechnen
            # Breiten vergleichen
            # if b_s < b_r : nach innen
            # if b_s > b_r : nach außen
            # else [0,0]
            polar_s = polarConverter2( s_area.p1, s_area.p2, snapshot.center)
            polar_r = polarConverter2( pair.p1, pair.p2, retina.center)
            breite_s = polar_s[1] - polar_s[0] 
            breite_r = polar_r[1] - polar_r[0]
            if breite_s < breite_r:
                translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)])*-1 )
            elif breite_s > breite_r:
                translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)]))
            else:
                translation_vectors.append(np.array([0,0]))

        # Berechnen der "opposite"-Vektoren mit denselben Methoden
        for s_area in snapshot.op_areas:
            abstand = []
            pair = 0
            #"Opposite"-Flächen des Snapshots werden nun mit den entsprechenden Retina-"Opposite"-Flächen zusammengepaart.
            for r_area in retina.op_areas:
                berechnung = abs(s_area.center - r_area.center)
                
                if (berechnung > (math.pi)):
                    berechnung = (2*math.pi - berechnung)
                abstand.append(abs(berechnung))
                
                if (abs(berechnung)) == min(abstand):
                    pair = r_area                         
            # Rotation Vektoren  
            if s_area.center < pair.center:
                rotation_vectors.append(rotateVector2DClockwise([math.cos(pair.center), math.sin(pair.center)]))
            elif s_area.center > pair.center:
                rotation_vectors.append(rotateVector2DCounterClockwise([math.cos(pair.center), math.sin(pair.center)]))
            else:
                rotation_vectors.append(np.array([0,0]))
            # Translation Vektoren
            # Breite der Areas berechnen
            # Breiten vergleichen
            # if b_s < b_r : nach innen
            # if b_s > b_r : nach außen
            # else [0,0]
            polar_s = polarConverter2( s_area.p1, s_area.p2, snapshot.center)
            polar_r = polarConverter2( pair.p1, pair.p2, retina.center)
            breite_s = polar_s[1] - polar_s[0] 
            breite_r = polar_r[1] - polar_r[0]
            if breite_s < breite_r:
                translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)])*-1 )
            elif breite_s > breite_r:
                translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)]))
            else:
                translation_vectors.append(np.array([0,0]))

        # Alle Rotations- und Translationsvektoren werden zusammengefügt.
        vt = np.array([0,0])
        vp = np.array([0,0])
        for x in rotation_vectors:
            vt = vt + x
        for x in translation_vectors:
            vp = vp + x

        '''
        Homing-Vektor ergibt sich aus vt + vp*3
        '''
        v = vt + (3 * vp) 
        
        nullpunkt_v = -(retina.center)
        thetaDiff = abs(polarConverter3(v, retina.center) - polarConverter3(nullpunkt_v, retina.center))
        avg_diff = np.rad2deg(thetaDiff)
        if (avg_diff > 180):
            avg_diff = 360 - avg_diff
        
        print("Abweichung in Grad: " + str(avg_diff))
        print("V = " + str(v)) #Ausgabe des Homing-Vektors
        
        '''
        Quiver-Plot wird erstellt
        '''
        ax.quiver(retina.center[0], retina.center[1], v[0], v[1])
        ax.set_title('ATS - Snapshotmodel')
        plt.show()

    elif eingabe == 2:
        #Für jede Koordinate von [-7;-7] bis [7;7] wird ein Homing-Vektor erstellt
        avg_diff_arr = [] # zur Bestimmung der Durchschnittsabweichung der Homing-Vektoren in Grad
        
        for j in range(15):
            for l in range(15):
                # Homing-Vektoren werden NICHT bei der Snapshot-Koordinate [0;0] oder bei den Landmarks erstellt.
                if ((j-7) == 0 and (l-7) == 0) or ((j-7) == LM1[0] and (l-7) == LM1[1]) or ((j-7) == LM2[0] and (l-7) == LM2[1]) or ((j-7) == LM3[0] and (l-7) == LM3[1]): 
                    continue
                
                '''Retina-Abbildung wird mit den entsprechenden Koordinaten erstellt.'''
                retina = takeSnapshot(np.array([j-7, l-7]), 2)
                
                '''
                Rotations- und Translations-Vektoren werden durch Zuordnung und Vergleich der Flächen bestimmt.
                '''
                rotation_vectors = []
                translation_vectors = []
                
                for s_area in snapshot.areas:
                    abstand = []
                    pair = 0
                    #Retina-Abbildungsflächen werden den passenden Snapshot-Flächen zugeordnet
                    for r_area in retina.areas:
                        berechnung = abs(s_area.center - r_area.center)
                        
                        if (berechnung > (math.pi)):
                            berechnung = (2*math.pi - berechnung)
                        abstand.append(abs(berechnung))
                        
                        if (abs(berechnung)) == min(abstand):
                            pair = r_area
                            
                    # Bestimmung der Rotationsvektoren
                    if s_area.center < pair.center:
                        #print(r_area.center)
                        rotation_vectors.append(rotateVector2DClockwise([math.cos(pair.center), math.sin(pair.center)]))
                    elif s_area.center > pair.center:
                        rotation_vectors.append(rotateVector2DCounterClockwise([math.cos(pair.center), math.sin(pair.center)]))
                    else:
                        rotation_vectors.append(np.array([0,0]))
                    # Determine Translation Vector
                    # Breite der Areas berechnen
                    # Breiten vergleichen
                    # if b_s < b_r : nach innen
                    # if b_s > b_r : nach außen
                    # else [0,0]
                    polar_s = polarConverter2( s_area.p1, s_area.p2, snapshot.center)
                    polar_r = polarConverter2( pair.p1, pair.p2, retina.center)
                    breite_s = polar_s[1] - polar_s[0] 
                    breite_r = polar_r[1] - polar_r[0]
                    if breite_s < breite_r:
                        translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)])*-1 )
                    elif breite_s > breite_r:
                        translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)]))
                    else:
                        translation_vectors.append(np.array([0,0]))

                # "Opposite"-Flächen des Snapshots werden nun mit den entsprechenden Retina-"Opposite"-Flächen zusammengepaart.
                for s_area in snapshot.op_areas:
                    abstand = []
                    pair = 0
                    # Find pair
                    for r_area in retina.op_areas:
                        berechnung = abs(s_area.center - r_area.center)
                        
                        if (berechnung > (math.pi)):
                            berechnung = (2*math.pi - berechnung)
                        abstand.append(abs(berechnung))
                        
                        if (abs(berechnung)) == min(abstand):
                            pair = r_area                  
                            
                    # Bestimmung der Rotationsvektoren
                    if s_area.center < pair.center:
                        rotation_vectors.append(rotateVector2DClockwise([math.cos(pair.center), math.sin(pair.center)]))
                    elif s_area.center > pair.center:
                        rotation_vectors.append(rotateVector2DCounterClockwise([math.cos(pair.center), math.sin(pair.center)]))
                    else:
                        rotation_vectors.append(np.array([0,0]))
                        
                    # Bestimmung der Translationsvektoren
                    polar_s = polarConverter2( s_area.p1, s_area.p2, snapshot.center)
                    polar_r = polarConverter2( pair.p1, pair.p2, retina.center)
                    breite_s = polar_s[1] - polar_s[0] 
                    breite_r = polar_r[1] - polar_r[0]
                    if breite_s < breite_r:
                        translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)])*-1 )
                    elif breite_s > breite_r:
                        translation_vectors.append(np.array([math.cos(pair.center), math.sin(pair.center)]))
                    else:
                        translation_vectors.append(np.array([0,0]))

                # Alle Rotations- und Translationsvektoren werden zusammengefügt.
                vt = np.array([0,0])
                vp = np.array([0,0])
                for x in rotation_vectors:
                    vt = vt + x
                for x in translation_vectors:
                    vp = vp + x

                '''
                Homing-Vektor ergibt sich aus vt + vp*3
                '''
                v = vt + (3 * vp)
                
                nullpunkt_v = -(retina.center)
                thetaDiff = abs(polarConverter3(v, retina.center) - polarConverter3(nullpunkt_v, retina.center))
                
                avg_diff = np.rad2deg(thetaDiff)
                if (avg_diff > 180):
                    avg_diff = 360 - avg_diff
                
                avg_diff_arr.append(avg_diff)
                
                #print(str(avg_diff) + " | Koordinaten (" + str(retina.center) + ")")
                
                '''
                Quiver-Plot wird erstellt
                '''
                ax.quiver(retina.center[0], retina.center[1], v[0], v[1])
                #Displaying the plot
                
        print(str(np.average(avg_diff_arr)) + " => Durchschnittliche Abweichung in Grad")
        ax.set_title('ATS - Snapshotmodel')
        plt.show()