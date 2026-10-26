"""Ask AREE: a public chat that answers plain-language questions from AREE's evidence.

Claude answers by calling the read-only query layer (`query`) as tools, so its
numbers come from AREE's outputs. `chat.server` is the web service, which holds
the API key and enforces rate and spend limits; `chat.conversation` runs one turn;
`chat.tools` adapts the query functions for a public audience. See docs/ask_aree.md.

Requires the optional `chat` extra (Python >= 3.10): pip install -e ".[chat]"
"""
