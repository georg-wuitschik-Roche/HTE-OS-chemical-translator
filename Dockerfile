# Use the official lightweight Python image.
FROM python:3.9-slim
# Allow statements and log 
ENV PYTHONUNBUFFERED True
# Copy local code to the container image.
ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

# Install production dependencies.
RUN apt update
RUN apt-get -y install gcc
RUN apt-get -y install g++
RUN apt-get install -yqq libsm6 libxext6 libxrender-dev
RUN pip install -r requirements.txt
# Run
CMD ["python", "main.py"]
