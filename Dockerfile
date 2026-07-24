# Use a Python slim image as the base
FROM python:3.12-slim
 
# Install Node.js and npm
RUN apt-get update && apt-get install -y \
    nodejs \
    npm \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
 
# Set working directory
WORKDIR /app
 
# Install Python dependencies
RUN pip install --no-cache-dir \
    markdown2 \
    pygments \
    beautifulsoup4
 
# Install katex via npm
RUN npm install katex
 
# Verify katex installation
RUN npx katex --version
 
# Copy the converter script (if present in build context)
# COPY MDLicious2.py /app/MDLicious2.py
 
# Config path is passed as an argument at runtime:
# docker run --rm -v /path/to/config.json:/config.json mdlicious2 /config.json
ENTRYPOINT ["python3", "/app/MDLicious2.py"]