#!/usr/local/bin/python
""" This function is used to send emails to a recepient.

This is used by me to send an email to myself indicating that the script 
I'm running has completed.
Example:
sendMessage( joe@gmail.com, hector@whocares.com, "You suck.", usr, !PasWrd?)
"""
import smtplib
import sys

def sendMessage(addr_from, addr_to, message, username, password):
	# Send the damn message
	server = smtplib.SMTP('smtp.gmail.com:587')
	server.starttls()
	server.login(username,password)
	server.sendmail(email, email, message)
	server.quit()

if __name__ == '__main__':
# Grab args
	if len(sys.argv) != 6:
		sys.exit("Usage: %s [Email from:] [Email to:] [Message] [Username] [Password]" % sys.argv[0])
	else:
		sendMessage(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
