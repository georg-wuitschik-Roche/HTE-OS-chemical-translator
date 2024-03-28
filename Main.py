import uvicorn, os
from translator import app

## RUN!!!
if __name__ == "__main__":
    uvicorn.run(app, port=int(os.environ.get("PORT", 8000)), host="0.0.0.0")

# if __name__ == "__main__":
#     import socket

#     hostname = socket.gethostname()
#     # host = hostname + ".kau.roche.com"
#     host = hostname
#     # host="0.0.0.0"
#     uvicorn.run(app, host=host, port=3020)

