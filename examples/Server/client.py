#!/usr/bin/env python
# encoding: utf-8
"""
client.py

A simple client to the PopulateReactionServer

Created by Richard West on 2011-06-21.
Copyright (c) 2011 MIT. All rights reserved.
"""
import socket
client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
client_socket.connect(("localhost", 5000))

request = """
C2H4 (molecule/cm3) 0.01 
1 C 0 {2,D}
2 C 0 {1,D}

O (molecule/cm3) 0.01
1 O 2T

END
"""

print "SENDING REQUEST"
client_socket.send(request)
print "REQUEST SENT"

partial_response = client_socket.recv(4096)
response = partial_response
print "RECEIVING RESPONSE"
while partial_response:
    partial_response = client_socket.recv(4096)
    response += partial_response

print "RECEIVED RESPONSE"
print response
print "CLOSING SOCKET"
client_socket.close()
