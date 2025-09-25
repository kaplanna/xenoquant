import torch
from torch import nn
from remora.activations import swish
from remora import constants

# Main network class without positional encoding
class network(nn.Module):
    _variable_width_possible = False

    def __init__(
        self,
        size=constants.DEFAULT_NN_SIZE,  # Hidden size (used for d_model)
        kmer_len=constants.DEFAULT_KMER_LEN,
        num_out=2,  # Number of output classes
    ):
        super().__init__()

        # Convolution layers for signal (sigs)
        self.sig_conv1 = nn.Conv1d(1, 4, 5)
        self.sig_bn1 = nn.BatchNorm1d(4)
        self.sig_conv2 = nn.Conv1d(4, 16, 5)
        self.sig_bn2 = nn.BatchNorm1d(16)
        self.sig_conv3 = nn.Conv1d(16, size, 9, 3)
        self.sig_bn3 = nn.BatchNorm1d(size)

        # Convolution layers for sequence (seqs)
        self.seq_conv1 = nn.Conv1d(kmer_len * 4, 16, 5)
        self.seq_bn1 = nn.BatchNorm1d(16)
        self.seq_conv2 = nn.Conv1d(16, size, 13, 3)
        self.seq_bn2 = nn.BatchNorm1d(size)

        # Merge conv layer for sigs and seqs
        self.merge_conv1 = nn.Conv1d(size * 2, size, 5)
        self.merge_bn = nn.BatchNorm1d(size)

        # Transformer encoder layers
        self.transformer_layer = nn.TransformerEncoderLayer(d_model=size, nhead=8)
        self.transformer_encoder = nn.TransformerEncoder(self.transformer_layer, num_layers=2)

        # Fully connected output layer
        self.fc = nn.Linear(size, num_out)

        # Dropout layer
        self.dropout = nn.Dropout(p=0.3)

    def forward(self, sigs, seqs):
        # Process signal (sigs) through convolutional layers
        sigs_x = swish(self.sig_bn1(self.sig_conv1(sigs)))
        sigs_x = swish(self.sig_bn2(self.sig_conv2(sigs_x)))
        sigs_x = swish(self.sig_bn3(self.sig_conv3(sigs_x)))

        # Process basecall sequences (seqs) through convolutional layers
        seqs_x = swish(self.seq_bn1(self.seq_conv1(seqs)))
        seqs_x = swish(self.seq_bn2(self.seq_conv2(seqs_x)))

        # Concatenate sigs and seqs along the feature dimension
        z = torch.cat((sigs_x, seqs_x), 1)

        # Process through merged conv layers
        z = swish(self.merge_bn(self.merge_conv1(z)))

        # Permute to (seq_len, batch_size, d_model) for transformer
        z = z.permute(2, 0, 1)

        # Pass through transformer encoder layers
        z = self.transformer_encoder(z)

        # Take the last time step (seq_len dimension is reduced to 1)
        z = z[-1]

        # Permute back to (batch_size, d_model)
        z = z.permute(0, 1)

        # Fully connected layer to produce final output
        z = self.fc(z)

        return z

