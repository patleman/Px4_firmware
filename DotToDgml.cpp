
//Main agenda of this file is to create a new dgml file from dot file 
/*
    example in dgml;

    for node:
    <Node Id="FlightTasks" Category="Class" Bounds="1803.70443333333,429.826452367872,76.7933333333333,25.96" File="D:\git\px4\Firmware\src\lib\FlightTasks\SubscriptionArray.hpp" />
    <Node/>
    attribute:
    1)Id 2)Category 3)Bounds 4)File

    for link:
    <Link Source="FlightTasks" Target="vehicle_command_ack" Bounds="1880.49776666667,443.997056009539,459.761590279912,14.2562850165528" />
    attribute :
    1) Source 2)Target 3)Bounds


    example in dot file

    module
    m_land_detector [label=land_detector color="#666666" fontcolor="#ffffff" fontsize=16 shape=box style=filled]
    topic
    t_position_setpoint_triplet [label=position_setpoint_triplet color="#5841d8" fontcolor="#ffffff" shape=ellipse style=filled]

        
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

using namespace std;

int main(){
    /*
    general approach : first read dot file and simultaneously 
 
    */
    // new file
    FILE *fwptr;
    fwptr=fopen("px4created.dgml","w");

    char wcontent[201222];

    // writing first line
    char first_line[]="<\?xml version=\"1.0\" encoding=\"utf-8\"\?>\n";
    strcpy(wcontent,first_line);
     

    
    
    fwrite ( wcontent, sizeof(char), sizeof(first_line), fwptr);
    fclose (fwptr);

    // reading file
   // FILE *frptr;
   // frptr=fopen("graph.fv","r");




}

