# Chemical Translator on Google Cloud Run

Adapted from Torsten Schindler's code, designed to run in a Google Cloud Run Container. The core of this code was built by Torsten Schindler. It was later adapted for use on the Google Cloud Platform and amended with password handling by Georg Wuitschik.

## Getting started

In order to prevent others from using your service, choose a password and replace all occurrances of _<\<API-Password\>>_ in Translator.py You will need this password to get Google Sheets and Spotfire to run.

Follow the steps outlined in this web article to get going: (https://medium.com/codex/secured-serverless-fastapi-with-google-cloud-run-66242b916b46)

After successfully setting up the API, you'll get the URL of the service which you need also for Google Sheets and Spotfire to run (placeholder there is <\<Link to API\>> ). Check Google_Cloud_Run_URL.png to see where to find it.

It's hard and complicated, if you haven't done it before. Don't give up. 

In Google Sheets, this happens whenever a molecular drawing is generated for presentations and wherever things like Inchi or Inchi-Key had to be generated. 
In Spotfire, the API is used for creating isotope patterns of the reaction components. 

