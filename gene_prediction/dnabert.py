import logging

import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

tokenizer = AutoTokenizer.from_pretrained(
    "zhihan1996/DNABERT-2-117M", trust_remote_code=True
)
model = AutoModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)

# Example DNA sequence
sequence = "ATGCGTACGTAGCTAGCTAGCTAGC"

# Tokenize using BPE
tokens = tokenizer(sequence, return_tensors="pt")

# Run inference
with torch.no_grad():
    outputs = model(**tokens)

# Print the outputs to understand its structure
print(outputs)

# Assuming the first element of the tuple is the logits
logits = outputs[0]
predictions = torch.argmax(logits, dim=-1)

logging.info(predictions)
