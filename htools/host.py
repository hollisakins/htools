import socket
hostname = socket.gethostname()

if hostname == 'cns-s-pmaa65432':
    hostname = 'patrick' # desktop

# elif hostname == '':
#     hostname = 'gerald' # laptop


elif hostname.startswith('cns'):
    raise Exception("Not sure what machine you're on.")

else:
    hostname = 'candide'