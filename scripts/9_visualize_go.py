import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read your Excel file
# Upload your file to Colab first, then use the correct path
df = pd.read_excel('x.xlsx', sheet_name='GO terms unmapped reads')

# Sort by count and take top 20
top_terms = df.nlargest(20, 'COUNT of GO_Terms')

# Horizontal bar plot for better readability of long names
plt.figure(figsize=(12, 10))
sns.barplot(data=top_terms, y='Name', x='COUNT of GO_Terms', orient='h')
plt.title('Top 20 GO Terms by Count')
plt.xlabel('Count')
plt.ylabel('GO Term Name')
plt.tight_layout()
plt.show()
